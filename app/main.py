import streamlit as st

from parser import load_gene_summary
from plots import build_volcano_plot

st.set_page_config(page_title="MAGeCK Dashboard", layout="wide")
st.title("MAGeCK Dashboard")

uploaded_file = st.file_uploader("Upload gene_summary.txt", type=["txt", "tsv"])

if uploaded_file is not None:
    df = load_gene_summary(uploaded_file)

    col1, col2 = st.columns(2)
    fdr_threshold = col1.number_input("FDR threshold", value=0.05, min_value=0.0, max_value=1.0, step=0.01)
    lfc_threshold = col2.number_input("LFC threshold", value=1.0, min_value=0.0, step=0.1)

    fig = build_volcano_plot(df, fdr_threshold=fdr_threshold, lfc_threshold=lfc_threshold)
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Gene summary table")
    st.dataframe(df, use_container_width=True)
else:
    st.info("Upload a MAGeCK gene_summary.txt file to get started.")