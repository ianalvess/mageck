import streamlit as st
import pandas as pd
from pathlib import Path

from parser import load_gene_summary, load_sgrna_summary, load_countsummary

from plots import (
    build_volcano_plot,
    build_read_distribution_plot,
    compute_gini_index,
    build_mapping_summary_plot,
    build_replicate_correlation_plot,
    figure_to_png_bytes,
    build_html_report,
    compare_gene_summaries,
    build_screen_comparison_plot,
    build_zerocounts_plot,
)

from pipeline import run_mageck_count, run_mageck_test

st.set_page_config(page_title="MAGeCK Dashboard", layout="wide")
st.title("MAGeCK Dashboard")

# Named figure variables, so the full-report builder at the bottom can
# find each one regardless of which tab last ran.
volcano_fig = None
qc_fig = None
mapping_fig = None
replicate_fig = None
comparison_fig = None

tab_pipeline, tab_volcano, tab_qc, tab_mapping, tab_replicates, tab_comparison = st.tabs(
    ["Run Pipeline", "Volcano Plot", "QC Report", "Mapping QC", "Replicate Correlation", "Screen Comparison"]
)

with tab_pipeline:
    st.subheader("1. Select input files")
    st.caption("Files must already be inside the data/input folder (mounted into the container).")

    data_input_dir = Path("data/input")
    available_files = sorted(str(p.relative_to(data_input_dir)) for p in data_input_dir.rglob("*") if p.is_file())

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

        st.subheader("2. Define groups")
        control_labels = st.multiselect("Control sample labels", options=sample_labels)
        treatment_labels = st.multiselect("Treatment sample labels", options=sample_labels)

        overlap = set(control_labels) & set(treatment_labels)

        st.subheader("3. Run")
        output_prefix = st.text_input("Output prefix", value="mageck_run")

        ready_to_run = (
            len(control_labels) > 0
            and len(treatment_labels) > 0
            and not overlap
            and output_prefix.strip() != ""
        )

        if overlap:
            st.error(f"These labels appear in both groups: {', '.join(overlap)}. A sample can't be both control and treatment.")
        elif not ready_to_run:
            st.info("Select the library file, all FASTQ files, and assign at least one control and one treatment label to continue.")

        if st.button("Run MAGeCK pipeline", disabled=not ready_to_run):
            work_dir = Path("data/pipeline_runs") / output_prefix
            work_dir.mkdir(parents=True, exist_ok=True)

            fastq_paths = [data_input_dir / f for f in fastq_choices]
            library_path = data_input_dir / library_choice

            with st.spinner("Running mageck count... this can take a while for large FASTQ files."):
                count_result = run_mageck_count(
                    fastq_paths=fastq_paths,
                    sample_labels=sample_labels,
                    library_path=library_path,
                    output_dir=work_dir,
                    output_prefix=output_prefix,
                )

            if not count_result.success:
                st.error("mageck count failed. See log below.")
                st.code(count_result.stderr or count_result.stdout)
            else:
                st.success("mageck count completed.")
                with st.expander("mageck count log"):
                    st.code(count_result.stdout + count_result.stderr)

                with st.spinner("Running mageck test..."):
                    test_result = run_mageck_test(
                        count_table_path=work_dir / f"{output_prefix}.count.txt",
                        treatment_labels=treatment_labels,
                        control_labels=control_labels,
                        output_dir=work_dir,
                        output_prefix=f"{output_prefix}_test",
                    )

                if not test_result.success:
                    st.error("mageck test failed. See log below.")
                    st.code(test_result.stderr or test_result.stdout)
                else:
                    st.session_state["pipeline_gene_summary"] = work_dir / f"{output_prefix}_test.gene_summary.txt"
                    st.session_state["pipeline_sgrna_summary"] = work_dir / f"{output_prefix}_test.sgrna_summary.txt"
                    st.session_state["pipeline_countsummary"] = work_dir / f"{output_prefix}.countsummary.txt"
                    st.session_state["pipeline_count_table"] = work_dir / f"{output_prefix}.count.txt"

                    st.success("mageck test completed! Results saved to: " + str(work_dir))
                    with st.expander("mageck test log"):
                        st.code(test_result.stdout + test_result.stderr)
with tab_volcano:
    source = st.radio(
        "Data source",
        ["Upload file", "Use last pipeline run"],
        key="volcano_source",
        horizontal=True,
    )

    if source == "Use last pipeline run":
        if "pipeline_gene_summary" not in st.session_state:
            st.warning("No pipeline run yet. Run the pipeline in the 'Run Pipeline' tab first, or switch to 'Upload file'.")
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
        st.plotly_chart(volcano_fig, use_container_width=True)

        st.download_button(
            "Download volcano plot (PNG)",
            data=figure_to_png_bytes(volcano_fig),
            file_name="volcano_plot.png",
            mime="image/png",
        )

        st.subheader("Gene summary table")
        st.dataframe(gene_df, use_container_width=True)
    elif source == "Upload file":
        st.info("Upload a MAGeCK gene_summary.txt file to get started.")

with tab_qc:
    source = st.radio(
        "Data source",
        ["Upload file", "Use last pipeline run"],
        key="qc_source",
        horizontal=True,
    )

    if source == "Use last pipeline run":
        if "pipeline_sgrna_summary" not in st.session_state:
            st.warning("No pipeline run yet. Run the pipeline in the 'Run Pipeline' tab first, or switch to 'Upload file'.")
            sgrna_df = None
        else:
            sgrna_df = load_sgrna_summary(st.session_state["pipeline_sgrna_summary"])
            st.caption(f"Using: {st.session_state['pipeline_sgrna_summary']}")
    else:
        uploaded_sgrna_file = st.file_uploader("Upload sgrna_summary.txt", type=["txt", "tsv"], key="sgrna_summary")
        sgrna_df = load_sgrna_summary(uploaded_sgrna_file) if uploaded_sgrna_file is not None else None

    if sgrna_df is not None:
        numeric_columns = sgrna_df.select_dtypes(include="number").columns.tolist()
        count_column = st.selectbox("Column to use as read count", numeric_columns)

        gini = compute_gini_index(sgrna_df[count_column])
        st.metric("Gini index", f"{gini:.3f}")

        qc_fig = build_read_distribution_plot(sgrna_df, count_column)
        st.plotly_chart(qc_fig, use_container_width=True)

        st.download_button(
            "Download read distribution plot (PNG)",
            data=figure_to_png_bytes(qc_fig),
            file_name="read_distribution.png",
            mime="image/png",
        )

        st.subheader("sgRNA summary table")
        st.dataframe(sgrna_df, use_container_width=True)
    elif source == "Upload file":
        st.info("Upload a MAGeCK sgrna_summary.txt file to get started.")

with tab_mapping:
    source = st.radio(
        "Data source",
        ["Upload file", "Use last pipeline run"],
        key="mapping_source",
        horizontal=True,
    )

    if source == "Use last pipeline run":
        if "pipeline_countsummary" not in st.session_state:
            st.warning("No pipeline run yet. Run the pipeline in the 'Run Pipeline' tab first, or switch to 'Upload file'.")
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
        st.plotly_chart(mapping_fig, use_container_width=True)

        st.download_button(
            "Download mapping summary plot (PNG)",
            data=figure_to_png_bytes(mapping_fig),
            file_name="mapping_summary.png",
            mime="image/png",
        )

        zerocounts_fig = build_zerocounts_plot(count_df)
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

with tab_replicates:
    source = st.radio(
        "Data source",
        ["Upload file", "Use last pipeline run"],
        key="replicates_source",
        horizontal=True,
    )

    if source == "Use last pipeline run":
        if "pipeline_count_table" not in st.session_state:
            st.warning("No pipeline run yet. Run the pipeline in the 'Run Pipeline' tab first, or switch to 'Upload file'.")
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
        col1, col2 = st.columns(2)
        x_column = col1.selectbox("Replicate 1", numeric_columns, index=0)
        y_column = col2.selectbox("Replicate 2", numeric_columns, index=min(1, len(numeric_columns) - 1))

        if x_column == y_column:
            st.warning("Select two different columns to compare.")
        else:
            replicate_fig = build_replicate_correlation_plot(counts_df, x_column, y_column)
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

with tab_comparison:
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

st.divider()
st.subheader("Full report")

available_figures = {
    "Volcano Plot": volcano_fig,
    "QC Report": qc_fig,
    "Mapping QC": mapping_fig,
    "Replicate Correlation": replicate_fig,
    "Screen Comparison": comparison_fig,
}
num_available = sum(fig is not None for fig in available_figures.values())

if num_available == 0:
    st.info("Load data in at least one tab above to enable the full report download.")
else:
    st.write(f"{num_available} chart(s) ready to include in the report.")
    report_html = build_html_report(available_figures)
    st.download_button(
        "Download full HTML report",
        data=report_html,
        file_name="mageck_report.html",
        mime="text/html",
    )