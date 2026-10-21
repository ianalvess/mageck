import streamlit as st
import pandas as pd
from pathlib import Path

from parser import load_gene_summary, load_countsummary, load_mle_gene_summary

from plots import (
    build_volcano_plot,
    build_read_distribution_plot,
    compute_gini_index,
    build_mapping_summary_plot,
    build_replicate_correlation_plot,
    figure_to_png_bytes,
    compare_gene_summaries,
    build_screen_comparison_plot,
    build_zerocounts_plot,
    build_correlation_heatmap,
    build_pca_plot,
    build_beta_comparison_plot,
    build_beta_score_plot,
    build_gc_content_plot,
    build_quality_profile_plot,
)
from fastq_qc import compute_gc_distribution, compute_quality_profile
from library_check import filter_observed, summarize_observed

from pipeline import run_mageck_count, run_mageck_test, run_mageck_mle

st.set_page_config(page_title="MAGeCK Dashboard", layout="wide")
st.title("MAGeCK Dashboard")

DATA_INPUT_DIR = Path("data/input")


def list_available_files():
    """List all files under data/input, as paths relative to that folder.

    Excludes hidden files (starting with a dot), such as .gitkeep,
    which exist only to keep the folder tracked in git.
    """
    return sorted(
        str(p.relative_to(DATA_INPUT_DIR))
        for p in DATA_INPUT_DIR.rglob("*")
        if p.is_file() and not p.name.startswith(".")
    )


section = st.sidebar.radio(
    "Section",
    ["Pipeline", "Quality Control", "Results", "Comparisons"],
)

if section == "Pipeline":
    st.header("Pipeline")
    sub_count, sub_test, sub_mle, sub_fastqqc = st.tabs(["Run Count", "Run Test", "Run MLE", "Raw FASTQ QC"])

    with sub_count:
        st.write(
            "`mageck count` aligns raw FASTQ reads to your sgRNA library and builds a table of read "
            "counts per sgRNA per sample, along with per-sample QC metrics (mapping rate, Gini index, "
            "zero-count sgRNAs). This is the first step of any screen analysis, required before "
            "running Test or MLE."
        )
        st.caption("Files must already be inside the data/input folder (mounted into the container).")

        available_files = list_available_files()

        if not available_files:
            st.warning("No files found in data/input. Copy your FASTQ and library files there first.")
        else:
            library_choice = st.selectbox("sgRNA library file (CSV)", available_files, key="library_choice")

            num_samples = st.number_input("Number of FASTQ samples", min_value=2, max_value=10, value=2, step=1)

            fastq_choices = []
            sample_labels = []
            for i in range(int(num_samples)):
                col1, col2 = st.columns([2, 1])
                fastq_choice = col1.selectbox(f"FASTQ file #{i + 1}", available_files, key=f"fastq_choice_{i}")
                label = col2.text_input(f"Label #{i + 1}", value=f"sample{i + 1}", key=f"label_{i}")
                fastq_choices.append(fastq_choice)
                sample_labels.append(label)

            count_output_prefix = st.text_input("Output prefix", value="mageck_run", key="count_output_prefix")

            ready_to_run_count = bool(fastq_choices) and bool(library_choice) and count_output_prefix.strip() != ""

            if st.button("Run mageck count", disabled=not ready_to_run_count):
                count_work_dir = Path("data/pipeline_runs") / count_output_prefix
                count_work_dir.mkdir(parents=True, exist_ok=True)

                fastq_paths = [DATA_INPUT_DIR / f for f in fastq_choices]
                library_path = DATA_INPUT_DIR / library_choice

                with st.spinner("Running mageck count... this can take a while for large FASTQ files."):
                    count_result = run_mageck_count(
                        fastq_paths=fastq_paths,
                        sample_labels=sample_labels,
                        library_path=library_path,
                        output_dir=count_work_dir,
                        output_prefix=count_output_prefix,
                    )

                if not count_result.success:
                    st.error("mageck count failed. See log below.")
                    st.code(count_result.stderr or count_result.stdout)
                else:
                    st.session_state["pipeline_count_table"] = count_work_dir / f"{count_output_prefix}.count.txt"
                    st.session_state["pipeline_countsummary"] = count_work_dir / f"{count_output_prefix}.countsummary.txt"
                    st.session_state["last_count_labels"] = sample_labels

                    st.success("mageck count completed! Count table is ready for Run Test, Run MLE, or Comparisons.")
                    with st.expander("mageck count log"):
                        st.code(count_result.stdout + count_result.stderr)

    with sub_test:
        st.write(
            "`mageck test` uses the RRA algorithm to compare a treatment group against a control "
            "group, ranking genes by consistent depletion or enrichment of their sgRNAs. It produces "
            "the gene- and sgRNA-level significance results used in the Test Results plot. Best for a "
            "simple two-group comparison (one control vs. one treatment condition)."
        )

        available_files_test = list_available_files()

        test_count_source = st.radio(
            "Count table source",
            ["Use last count run", "Select from data/input"],
            key="test_count_source",
            horizontal=True,
        )

        if test_count_source == "Use last count run":
            if "pipeline_count_table" not in st.session_state:
                st.warning("No count table yet. Run mageck count in the Run Count tab first, or switch to 'Select from data/input'.")
                test_count_path = None
            else:
                test_count_path = st.session_state["pipeline_count_table"]
                st.caption(f"Using: {test_count_path}")
        else:
            test_count_choice = st.selectbox("Count table file", available_files_test, key="test_count_choice")
            test_count_path = DATA_INPUT_DIR / test_count_choice if test_count_choice else None

        available_labels = st.session_state.get("last_count_labels", [])
        if available_labels:
            control_labels = st.multiselect("Control sample labels", options=available_labels, key="test_control_labels")
            treatment_labels = st.multiselect("Treatment sample labels", options=available_labels, key="test_treatment_labels")
        else:
            st.info("Sample labels for this count table are unknown here. Enter them manually, comma-separated.")
            control_text = st.text_input("Control sample labels (comma-separated)", key="test_control_text")
            treatment_text = st.text_input("Treatment sample labels (comma-separated)", key="test_treatment_text")
            control_labels = [s.strip() for s in control_text.split(",") if s.strip()]
            treatment_labels = [s.strip() for s in treatment_text.split(",") if s.strip()]

        overlap = set(control_labels) & set(treatment_labels)
        test_output_prefix = st.text_input("Output prefix", value="mageck_test", key="test_output_prefix")

        ready_to_run_test = (
            test_count_path is not None
            and len(control_labels) > 0
            and len(treatment_labels) > 0
            and not overlap
            and test_output_prefix.strip() != ""
        )

        if overlap:
            st.error(f"These labels appear in both groups: {', '.join(overlap)}. A sample can't be both control and treatment.")
        elif not ready_to_run_test:
            st.info("Provide a count table and assign at least one control and one treatment label to continue.")

        if st.button("Run mageck test", disabled=not ready_to_run_test):
            test_work_dir = Path("data/pipeline_runs") / test_output_prefix
            test_work_dir.mkdir(parents=True, exist_ok=True)

            with st.spinner("Running mageck test..."):
                test_result = run_mageck_test(
                    count_table_path=test_count_path,
                    treatment_labels=treatment_labels,
                    control_labels=control_labels,
                    output_dir=test_work_dir,
                    output_prefix=test_output_prefix,
                )

            if not test_result.success:
                st.error("mageck test failed. See log below.")
                st.code(test_result.stderr or test_result.stdout)
            else:
                st.session_state["pipeline_gene_summary"] = test_work_dir / f"{test_output_prefix}.gene_summary.txt"
                st.session_state["pipeline_sgrna_summary"] = test_work_dir / f"{test_output_prefix}.sgrna_summary.txt"

                st.success("mageck test completed! Results ready in Results > Test Results.")
                with st.expander("mageck test log"):
                    st.code(test_result.stdout + test_result.stderr)

    with sub_mle:
        st.write(
            "`mageck mle` fits a statistical model across all samples at once, using a design matrix "
            "to describe your experimental conditions. It estimates a beta score per gene per "
            "condition, similar in spirit to the LFC from `mageck test`. Use this instead of "
            "`mageck test` when you have more than two groups, multiple treatments, or covariates "
            "you want to control for."
        )

        available_files_mle = list_available_files()

        if "pipeline_count_table" not in st.session_state and not available_files_mle:
            st.warning("No count table available. Run 'mageck count' in the Run Count tab first, or add files to data/input.")
        else:
            count_source = st.radio(
                "Count table source",
                ["Use last count run", "Select from data/input"],
                key="mle_count_source",
                horizontal=True,
            )

            if count_source == "Use last count run":
                if "pipeline_count_table" not in st.session_state:
                    st.warning("No count table yet.")
                    mle_count_path = None
                else:
                    mle_count_path = st.session_state["pipeline_count_table"]
                    st.caption(f"Using: {mle_count_path}")
            else:
                mle_count_choice = st.selectbox("Count table file", available_files_mle, key="mle_count_choice")
                mle_count_path = DATA_INPUT_DIR / mle_count_choice if mle_count_choice else None

            design_choice = st.selectbox("Design matrix file", available_files_mle, key="mle_design_choice")
            mle_output_prefix = st.text_input("Output prefix", value="mle_run", key="mle_output_prefix")

            mle_ready = mle_count_path is not None and design_choice is not None and mle_output_prefix.strip() != ""

            if st.button("Run MAGeCK MLE", disabled=not mle_ready):
                mle_work_dir = Path("data/pipeline_runs") / mle_output_prefix
                mle_work_dir.mkdir(parents=True, exist_ok=True)
                design_path = DATA_INPUT_DIR / design_choice

                with st.spinner("Running mageck mle..."):
                    mle_result = run_mageck_mle(
                        count_table_path=mle_count_path,
                        design_matrix_path=design_path,
                        output_dir=mle_work_dir,
                        output_prefix=mle_output_prefix,
                    )

                if not mle_result.success:
                    st.error("mageck mle failed. See log below.")
                    st.code(mle_result.stderr or mle_result.stdout)
                else:
                    st.session_state["pipeline_mle_gene_summary"] = mle_work_dir / f"{mle_output_prefix}.gene_summary.txt"
                    st.success("mageck mle completed! Go to Results > MLE Results to see the beta scores.")
                    with st.expander("mageck mle log"):
                        st.code(mle_result.stdout + mle_result.stderr)

    with sub_fastqqc:
        st.write(
            "Computes two standard sequencing QC metrics directly from a raw FASTQ file: the GC "
            "content distribution across reads, and the median base quality by position along the "
            "read. Use this to catch sequencing problems (e.g. low quality, unexpected GC skew) "
            "before running the rest of the pipeline."
        )

        available_files_qc = list_available_files()

        if not available_files_qc:
            st.warning("No files found in data/input. Copy your FASTQ files there first.")
        else:
            fastq_choice_qc = st.selectbox("FASTQ file", available_files_qc, key="fastqqc_choice")
            sample_name_qc = st.text_input("Sample name (for plot titles)", value=Path(fastq_choice_qc).stem)

            if st.button("Compute QC"):
                fastq_path_qc = DATA_INPUT_DIR / fastq_choice_qc

                with st.spinner(f"Reading {fastq_choice_qc} and computing GC content..."):
                    gc_df = compute_gc_distribution(fastq_path_qc)

                with st.spinner("Computing base quality profile..."):
                    quality_df = compute_quality_profile(fastq_path_qc)

                st.session_state["fastqqc_gc_df"] = gc_df
                st.session_state["fastqqc_quality_df"] = quality_df
                st.session_state["fastqqc_sample_name"] = sample_name_qc
                st.success(f"Processed {len(gc_df):,} reads from {fastq_choice_qc}.")

            if "fastqqc_gc_df" in st.session_state:
                fastq_gc_fig = build_gc_content_plot(
                    st.session_state["fastqqc_gc_df"], sample_name=st.session_state["fastqqc_sample_name"]
                )
                st.session_state["fig_fastq_gc"] = fastq_gc_fig
                st.plotly_chart(fastq_gc_fig, use_container_width=True)
                st.download_button(
                    "Download GC content plot (PNG)",
                    data=figure_to_png_bytes(fastq_gc_fig),
                    file_name="gc_content.png",
                    mime="image/png",
                )

                fastq_quality_fig = build_quality_profile_plot(
                    st.session_state["fastqqc_quality_df"], sample_name=st.session_state["fastqqc_sample_name"]
                )
                st.session_state["fig_fastq_quality"] = fastq_quality_fig
                st.plotly_chart(fastq_quality_fig, use_container_width=True)
                st.download_button(
                    "Download base quality plot (PNG)",
                    data=figure_to_png_bytes(fastq_quality_fig),
                    file_name="base_quality.png",
                    mime="image/png",
                )

elif section == "Quality Control":
    st.header("Quality Control")
    sub_mapping, sub_qc, sub_overview = st.tabs(["Mapping QC", "QC Report", "Sample Overview"])

    with sub_mapping:
        st.write(
            "Shows the read mapping rate and the number of sgRNAs with zero reads per sample, "
            "computed by `mageck count`. A low mapping rate or a high number of zero-count sgRNAs "
            "for a sample usually points to a sequencing, PCR, or library preparation problem."
        )

        source = st.radio(
            "Data source",
            ["Upload file", "Use last count run"],
            key="mapping_source",
            horizontal=True,
        )

        if source == "Use last count run":
            if "pipeline_countsummary" not in st.session_state:
                st.warning("No count table yet. Run mageck count in Pipeline > Run Count first, or switch to 'Upload file'.")
                count_df = None
            else:
                count_df = load_countsummary(st.session_state["pipeline_countsummary"])
                st.caption(f"Using: {st.session_state['pipeline_countsummary']}")
        else:
            uploaded_count_file = st.file_uploader("Upload countsummary.txt", type=["txt", "tsv"], key="countsummary")
            count_df = load_countsummary(uploaded_count_file) if uploaded_count_file is not None else None

        if count_df is not None:
            low_mapping = count_df[count_df["Percentage"] < 0.8]
            if not low_mapping.empty:
                st.warning(f"Samples below 80% mapping: {', '.join(low_mapping['Label'])}")

            mapping_fig = build_mapping_summary_plot(count_df)
            st.session_state["fig_mapping"] = mapping_fig
            st.plotly_chart(mapping_fig, use_container_width=True)
            st.download_button(
                "Download mapping summary plot (PNG)",
                data=figure_to_png_bytes(mapping_fig),
                file_name="mapping_summary.png",
                mime="image/png",
            )

            zerocounts_fig = build_zerocounts_plot(count_df)
            st.session_state["fig_zerocounts"] = zerocounts_fig
            st.plotly_chart(zerocounts_fig, use_container_width=True)
            st.download_button(
                "Download zero-count plot (PNG)",
                data=figure_to_png_bytes(zerocounts_fig),
                file_name="zerocounts.png",
                mime="image/png",
            )

            st.subheader("Count summary table")
            st.dataframe(count_df, use_container_width=True)
        elif source == "Upload file":
            st.info("Upload a MAGeCK countsummary.txt file to get started.")

    with sub_qc:
        st.write(
            "Computes the Gini index and read-count distribution for one individual sample at a "
            "time, from the raw per-sgRNA count table. A high Gini index (uneven read distribution) "
            "or a distribution with many low-count sgRNAs can signal library or sequencing problems "
            "for that specific sample."
        )

        source = st.radio(
            "Data source",
            ["Upload file", "Use last count run"],
            key="qc_source",
            horizontal=True,
        )

        if source == "Use last count run":
            if "pipeline_count_table" not in st.session_state:
                st.warning(
                    "No count table yet. Run mageck count in Pipeline > Run Count first, or switch to 'Upload file'.")
                counts_qc_df = None
            else:
                counts_qc_df = pd.read_csv(st.session_state["pipeline_count_table"], sep="\t")
                st.caption(f"Using: {st.session_state['pipeline_count_table']}")
        else:
            uploaded_counts_qc_file = st.file_uploader(
                "Upload a per-sample count table (e.g. count.txt)", type=["txt", "tsv"], key="qc_counts"
            )
            counts_qc_df = pd.read_csv(uploaded_counts_qc_file,
                                       sep="\t") if uploaded_counts_qc_file is not None else None

        if counts_qc_df is not None:
            sample_columns = counts_qc_df.select_dtypes(include="number").columns.tolist()
            id_column = "sgRNA" if "sgRNA" in counts_qc_df.columns else counts_qc_df.columns[0]

            observed_summary = summarize_observed(counts_qc_df, sample_columns, id_column=id_column)
            if observed_summary["unobserved"] > 0:
                st.warning(
                    f"{observed_summary['unobserved']:,} of {observed_summary['total']:,} sgRNAs "
                    f"({observed_summary['unobserved_fraction']:.1%}) have zero reads in every sample. "
                    "This usually means part of the library file (for example one half of a split "
                    "library) was not used in this screen."
                )
                with st.expander("Never-observed sgRNAs by ID prefix"):
                    st.dataframe(observed_summary["by_prefix"], use_container_width=True)

                if st.checkbox(
                        "Exclude never-observed sgRNAs from the analysis below",
                        value=True,
                        key="qc_exclude_unobserved",
                ):
                    counts_qc_df = filter_observed(counts_qc_df, sample_columns)
                    st.caption(f"Showing {len(counts_qc_df):,} sgRNAs after excluding never-observed ones.")

            count_column = st.selectbox("Sample", sample_columns)
            sample_counts = counts_qc_df[count_column]

            col1, col2, col3 = st.columns(3)
            col1.metric(
                "Gini (log2, as MAGeCK)",
                f"{compute_gini_index(sample_counts, log_scale=True):.3f}",
                help="Same definition as the GiniIndex in the MAGeCK countsummary (computed on log2(count + 1)).",
            )
            col2.metric(
                "Gini (raw counts)",
                f"{compute_gini_index(sample_counts):.3f}",
                help="Sensitive to a few dominant sgRNAs: values near 1 mean almost all reads sit in very few sgRNAs.",
            )
            col3.metric("sgRNAs with zero reads", f"{(sample_counts == 0).mean():.1%}")

            qc_fig = build_read_distribution_plot(counts_qc_df, count_column)
            st.session_state["fig_qc"] = qc_fig
            st.plotly_chart(qc_fig, use_container_width=True)
            st.download_button(
                "Download read distribution plot (PNG)",
                data=figure_to_png_bytes(qc_fig),
                file_name="read_distribution.png",
                mime="image/png",
            )

            st.subheader("Count table")
            st.dataframe(counts_qc_df, use_container_width=True)
        elif source == "Upload file":
            st.info("Upload a per-sample sgRNA count table to see the read distribution and Gini index.")

    with sub_overview:
        st.write(
            "Shows how all samples relate to each other at once: a correlation heatmap and a PCA "
            "plot. Useful to spot batch effects or mislabeled samples — replicates of the same "
            "condition should cluster together. Like the Replicate Correlation tab, this is a "
            "QC/reproducibility check, not biological validation of the screen."
        )

        uploaded_overview_file = st.file_uploader(
            "Upload a per-sample count table (e.g. count.txt)", type=["txt", "tsv"], key="overview_counts"
        )

        if uploaded_overview_file is not None:
            overview_df = pd.read_csv(uploaded_overview_file, sep="\t")
            numeric_columns = overview_df.select_dtypes(include="number").columns.tolist()

            selected_samples = st.multiselect(
                "Samples to include", numeric_columns, default=numeric_columns
            )
            log_transform_overview = st.checkbox(
                "Apply log2(count + 1) transform", value=True, key="overview_log_transform"
            )

            if len(selected_samples) < 2:
                st.warning("Select at least two samples to compare.")
            else:
                st.subheader("Correlation heatmap")
                overview_heatmap_fig = build_correlation_heatmap(
                    overview_df, selected_samples, log_transform=log_transform_overview
                )
                st.session_state["fig_overview_heatmap"] = overview_heatmap_fig
                st.plotly_chart(overview_heatmap_fig, use_container_width=True)
                st.download_button(
                    "Download correlation heatmap (PNG)",
                    data=figure_to_png_bytes(overview_heatmap_fig),
                    file_name="correlation_heatmap.png",
                    mime="image/png",
                )

                st.subheader("PCA")
                overview_pca_fig = build_pca_plot(
                    overview_df, selected_samples, log_transform=log_transform_overview
                )
                st.session_state["fig_overview_pca"] = overview_pca_fig
                st.plotly_chart(overview_pca_fig, use_container_width=True)
                st.download_button(
                    "Download PCA plot (PNG)",
                    data=figure_to_png_bytes(overview_pca_fig),
                    file_name="pca_plot.png",
                    mime="image/png",
                )
        else:
            st.info("Upload a per-sample sgRNA count table to see the sample overview.")

elif section == "Results":
    st.header("Results")
    sub_volcano, sub_mle_results = st.tabs(["Test Results", "MLE Results"])

    with sub_volcano:
        source = st.radio(
            "Data source",
            ["Upload file", "Use last test run"],
            key="volcano_source",
            horizontal=True,
        )

        if source == "Use last test run":
            if "pipeline_gene_summary" not in st.session_state:
                st.warning("No test run yet. Run mageck test in Pipeline > Run Test first, or switch to 'Upload file'.")
                gene_df = None
            else:
                gene_df = load_gene_summary(st.session_state["pipeline_gene_summary"])
                st.caption(f"Using: {st.session_state['pipeline_gene_summary']}")
        else:
            uploaded_gene_file = st.file_uploader("Upload gene_summary.txt", type=["txt", "tsv"], key="gene_summary")
            gene_df = load_gene_summary(uploaded_gene_file) if uploaded_gene_file is not None else None

        if gene_df is not None:
            col1, col2 = st.columns(2)
            fdr_threshold = col1.number_input("FDR threshold", value=0.05, min_value=0.0, max_value=1.0, step=0.01)
            lfc_threshold = col2.number_input("LFC threshold", value=1.0, min_value=0.0, step=0.1)

            volcano_fig = build_volcano_plot(gene_df, fdr_threshold=fdr_threshold, lfc_threshold=lfc_threshold)
            st.session_state["fig_volcano"] = volcano_fig
            st.plotly_chart(volcano_fig, use_container_width=True)
            st.download_button(
                "Download test results plot (PNG)",
                data=figure_to_png_bytes(volcano_fig),
                file_name="test_results_plot.png",
                mime="image/png",
            )

            st.subheader("Gene summary table")
            st.dataframe(gene_df, use_container_width=True)
        elif source == "Upload file":
            st.info("Upload a MAGeCK gene_summary.txt file to get started.")

    with sub_mle_results:
        mle_source = st.radio(
            "Data source",
            ["Upload file", "Use last MLE run"],
            key="mle_view_source",
            horizontal=True,
        )

        if mle_source == "Use last MLE run":
            if "pipeline_mle_gene_summary" not in st.session_state:
                st.warning("No MLE run yet. Run MLE in Pipeline > Run MLE first, or switch to 'Upload file'.")
                mle_df, mle_conditions = None, []
            else:
                mle_df, mle_conditions = load_mle_gene_summary(st.session_state["pipeline_mle_gene_summary"])
                st.caption(f"Using: {st.session_state['pipeline_mle_gene_summary']}")
        else:
            uploaded_mle_file = st.file_uploader(
                "Upload MLE gene_summary.txt", type=["txt", "tsv"], key="mle_upload"
            )
            if uploaded_mle_file is not None:
                mle_df, mle_conditions = load_mle_gene_summary(uploaded_mle_file)
            else:
                mle_df, mle_conditions = None, []

        if mle_df is not None:
            fdr_methods = {"Permutation (fdr)": "fdr"}
            if any(col.endswith("|wald-fdr") for col in mle_df.columns):
                fdr_methods["Wald (wald-fdr)"] = "wald-fdr"

            fdr_method_label = st.radio("FDR method", list(fdr_methods), horizontal=True, key="mle_fdr_method")
            fdr_suffix = fdr_methods[fdr_method_label]
            if fdr_suffix == "wald-fdr":
                st.warning(
                    "Wald FDR can be over-optimistic with few replicates or sparse counts. "
                    "Confirm hits with another method (e.g. mageck test)."
                )

            mle_fdr_threshold = st.number_input(
                "FDR threshold", value=0.05, min_value=0.0, max_value=1.0, step=0.01, key="mle_fdr"
            )

            def count_below_threshold(condition):
                fdr_values = mle_df[f"{condition}|{fdr_suffix}"]
                return int((fdr_values < mle_fdr_threshold).sum()), int(fdr_values.notna().sum())

            if len(mle_conditions) < 2:
                single_condition = mle_conditions[0]
                n_sig, n_total = count_below_threshold(single_condition)
                st.metric(f"Genes with FDR < {mle_fdr_threshold} ({single_condition})", f"{n_sig:,} of {n_total:,}")
                st.caption(
                    f"This result has a single condition ({single_condition}), so beta scores are plotted "
                    "against significance. With two or more conditions, a beta vs. beta comparison is shown."
                )
                mle_fig = build_beta_score_plot(
                    mle_df, condition=single_condition, fdr_threshold=mle_fdr_threshold, fdr_suffix=fdr_suffix
                )
                st.session_state["fig_mle"] = mle_fig
                st.plotly_chart(mle_fig, use_container_width=True)
                st.download_button(
                    "Download beta score plot (PNG)",
                    data=figure_to_png_bytes(mle_fig),
                    file_name=f"beta_score_{single_condition}.png",
                    mime="image/png",
                )
            else:
                col1, col2 = st.columns(2)
                condition_x = col1.selectbox("X axis condition", mle_conditions, index=0)
                condition_y = col2.selectbox(
                    "Y axis condition", mle_conditions, index=min(1, len(mle_conditions) - 1)
                )

                if condition_x == condition_y:
                    st.warning("Select two different conditions to compare.")
                else:
                    metric_col1, metric_col2 = st.columns(2)
                    for metric_col, condition in ((metric_col1, condition_x), (metric_col2, condition_y)):
                        n_sig, n_total = count_below_threshold(condition)
                        metric_col.metric(f"Genes with FDR < {mle_fdr_threshold} ({condition})", f"{n_sig:,} of {n_total:,}")

                    mle_fig = build_beta_comparison_plot(
                        mle_df,
                        condition_x=condition_x,
                        condition_y=condition_y,
                        fdr_threshold=mle_fdr_threshold,
                        fdr_suffix=fdr_suffix,
                    )
                    st.session_state["fig_mle"] = mle_fig
                    st.plotly_chart(mle_fig, use_container_width=True)
                    st.download_button(
                        "Download beta comparison plot (PNG)",
                        data=figure_to_png_bytes(mle_fig),
                        file_name=f"beta_comparison_{condition_x}_vs_{condition_y}.png",
                        mime="image/png",
                    )

            st.subheader("MLE gene summary table")
            st.dataframe(mle_df, use_container_width=True)
        elif mle_source == "Upload file":
            st.info("Upload a MAGeCK MLE gene_summary.txt file to get started.")

elif section == "Comparisons":
    st.header("Comparisons")
    sub_replicates, sub_screen_comparison = st.tabs(["Replicate Correlation", "Screen Comparison"])

    with sub_replicates:
        st.write(
            "Compares raw sgRNA read counts between two samples to check reproducibility. This is "
            "a QC/reproducibility metric, not evidence on its own that the screen worked "
            "biologically — low correlation between replicates of the same condition usually "
            "signals a technical problem (pipetting, contamination, low coverage)."
        )

        source = st.radio(
            "Data source",
            ["Upload file", "Use last count run"],
            key="replicates_source",
            horizontal=True,
        )

        if source == "Use last count run":
            if "pipeline_count_table" not in st.session_state:
                st.warning("No count table yet. Run mageck count in Pipeline > Run Count first, or switch to 'Upload file'.")
                counts_df = None
            else:
                counts_df = pd.read_csv(st.session_state["pipeline_count_table"], sep="\t")
                st.caption(f"Using: {st.session_state['pipeline_count_table']}")
        else:
            uploaded_counts_file = st.file_uploader(
                "Upload a per-sample count table (e.g. count.txt)", type=["txt", "tsv"], key="counts"
            )
            counts_df = pd.read_csv(uploaded_counts_file, sep="\t") if uploaded_counts_file is not None else None

        if counts_df is not None:
            numeric_columns = counts_df.select_dtypes(include="number").columns.tolist()

            group_text = st.text_input(
                "Group each sample, comma-separated in the same order as the columns below "
                f"({', '.join(numeric_columns)}). Leave blank to skip grouping.",
                key="replicate_groups",
            )
            groups = [g.strip() for g in group_text.split(",")] if group_text.strip() else []

            suggest_pairs = len(groups) == len(numeric_columns) and len(groups) > 0
            if group_text.strip() and not suggest_pairs:
                st.warning(f"Expected {len(numeric_columns)} group labels, got {len(groups)}. Ignoring grouping.")

            if suggest_pairs:
                column_to_group = dict(zip(numeric_columns, groups))
                same_group_pairs = [
                    (a, b)
                    for i, a in enumerate(numeric_columns)
                    for b in numeric_columns[i + 1 :]
                    if column_to_group[a] == column_to_group[b]
                ]
                if same_group_pairs:
                    pair_labels = [f"{a} vs {b}" for a, b in same_group_pairs]
                    chosen_pair = st.selectbox("Suggested same-group pairs", pair_labels)
                    x_column, y_column = same_group_pairs[pair_labels.index(chosen_pair)]
                else:
                    st.info("No two samples share the same group label; falling back to manual selection.")
                    col1, col2 = st.columns(2)
                    x_column = col1.selectbox("Replicate 1", numeric_columns, index=0)
                    y_column = col2.selectbox("Replicate 2", numeric_columns, index=min(1, len(numeric_columns) - 1))
            else:
                col1, col2 = st.columns(2)
                x_column = col1.selectbox("Replicate 1", numeric_columns, index=0)
                y_column = col2.selectbox("Replicate 2", numeric_columns, index=min(1, len(numeric_columns) - 1))

            log_transform = st.checkbox("Apply log2(count + 1) transform", value=True)

            if x_column == y_column:
                st.warning("Select two different columns to compare.")
            else:
                replicate_fig = build_replicate_correlation_plot(
                    counts_df, x_column, y_column, log_transform=log_transform
                )
                st.session_state["fig_replicate"] = replicate_fig
                st.plotly_chart(replicate_fig, use_container_width=True)
                st.download_button(
                    "Download replicate correlation plot (PNG)",
                    data=figure_to_png_bytes(replicate_fig),
                    file_name="replicate_correlation.png",
                    mime="image/png",
                )

            st.subheader("Count table")
            st.dataframe(counts_df, use_container_width=True)
        elif source == "Upload file":
            st.info("Upload a per-sample sgRNA count table to compare replicates.")

    with sub_screen_comparison:
        st.write("Upload two gene_summary.txt files from different screens to compare gene-level results.")

        col1, col2 = st.columns(2)
        file1 = col1.file_uploader("Screen 1: gene_summary.txt", type=["txt", "tsv"], key="screen1")
        label1 = col1.text_input("Screen 1 label", value="Screen 1")
        file2 = col2.file_uploader("Screen 2: gene_summary.txt", type=["txt", "tsv"], key="screen2")
        label2 = col2.text_input("Screen 2 label", value="Screen 2")

        if file1 is not None and file2 is not None:
            comparison_df1 = load_gene_summary(file1)
            comparison_df2 = load_gene_summary(file2)

            fdr_threshold_comp = st.number_input(
                "FDR threshold for significance", value=0.05, min_value=0.0, max_value=1.0, step=0.01, key="comp_fdr"
            )

            merged = compare_gene_summaries(
                comparison_df1, comparison_df2, fdr_threshold=fdr_threshold_comp, label1=label1, label2=label2
            )

            st.subheader("Concordance summary")
            st.write(merged["concordance"].value_counts())

            comparison_fig = build_screen_comparison_plot(merged, label1=label1, label2=label2)
            st.session_state["fig_comparison"] = comparison_fig
            st.plotly_chart(comparison_fig, use_container_width=True)
            st.download_button(
                "Download comparison plot (PNG)",
                data=figure_to_png_bytes(comparison_fig),
                file_name="screen_comparison.png",
                mime="image/png",
            )

            st.subheader("Merged gene table")
            st.dataframe(merged, use_container_width=True)
        else:
            st.info("Upload both gene_summary.txt files to compare.")