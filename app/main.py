import streamlit as st

from parser import load_gene_summary, load_sgrna_summary, load_countsummary
from plots import (
    build_volcano_plot,
    build_read_distribution_plot,
    compute_gini_index,
    build_mapping_summary_plot,
)

st.set_page_config(page_title="MAGeCK Dashboard", layout="wide")
st.title("MAGeCK Dashboard")

tab_volcano, tab_qc, tab_mapping = st.tabs(["Volcano Plot", "QC Report", "Mapping QC"])

with tab_volcano:
    uploaded_gene_file = st.file_uploader("Upload gene_summary.txt", type=["txt", "tsv"], key="gene_summary")

    if uploaded_gene_file is not None:
        gene_df = load_gene_summary(uploaded_gene_file)

        col1, col2 = st.columns(2)
        fdr_threshold = col1.number_input("FDR threshold", value=0.05, min_value=0.0, max_value=1.0, step=0.01)
        lfc_threshold = col2.number_input("LFC threshold", value=1.0, min_value=0.0, step=0.1)

        fig = build_volcano_plot(gene_df, fdr_threshold=fdr_threshold, lfc_threshold=lfc_threshold)
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Gene summary table")
        st.dataframe(gene_df, use_container_width=True)
    else:
        st.info("Upload a MAGeCK gene_summary.txt file to get started.")

with tab_qc:
    uploaded_sgrna_file = st.file_uploader("Upload sgrna_summary.txt", type=["txt", "tsv"], key="sgrna_summary")

    if uploaded_sgrna_file is not None:
        sgrna_df = load_sgrna_summary(uploaded_sgrna_file)

        numeric_columns = sgrna_df.select_dtypes(include="number").columns.tolist()
        count_column = st.selectbox("Column to use as read count", numeric_columns)

        gini = compute_gini_index(sgrna_df[count_column])
        st.metric("Gini index", f"{gini:.3f}")

        fig = build_read_distribution_plot(sgrna_df, count_column)
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("sgRNA summary table")
        st.dataframe(sgrna_df, use_container_width=True)
    else:
        st.info("Upload a MAGeCK sgrna_summary.txt file to get started.")

with tab_mapping:
    uploaded_count_file = st.file_uploader("Upload countsummary.txt", type=["txt", "tsv"], key="countsummary")

    if uploaded_count_file is not None:
        count_df = load_countsummary(uploaded_count_file)

        low_mapping = count_df[count_df["Percentage"] < 0.8]
        if not low_mapping.empty:
            st.warning(f"Samples below 80% mapping: {', '.join(low_mapping['Label'])}")

        fig = build_mapping_summary_plot(count_df)
        st.plotly_chart(fig, use_container_width=True)

        st.subheader("Count summary table")
        st.dataframe(count_df, use_container_width=True)
    else:
        st.info("Upload a MAGeCK countsummary.txt file to get started.")