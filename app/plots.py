"""
Plotting functions for the MAGeCK dashboard.

These functions build Plotly figures from the DataFrames produced by
app.parser. Kept separate from the Streamlit app so they can be
tested and reused independently.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def build_volcano_plot(
    df: pd.DataFrame,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
) -> go.Figure:
    """Build an interactive volcano plot from a gene_summary DataFrame.

    Args:
        df: DataFrame produced by app.parser.load_gene_summary
            (must contain "id", "best_lfc", "best_fdr", "direction").
        fdr_threshold: FDR cutoff used to color genes as significant.
        lfc_threshold: absolute LFC cutoff used to color genes as significant.

    Returns:
        A Plotly Figure ready to be rendered (e.g. with st.plotly_chart).
    """
    plot_df = df.copy()

    plot_df["neg_log10_fdr"] = -np.log10(plot_df["best_fdr"].clip(lower=1e-300))

    is_significant = (plot_df["best_fdr"] < fdr_threshold) & (
        plot_df["best_lfc"].abs() > lfc_threshold
    )
    plot_df["significance"] = is_significant.map({True: "Significant", False: "Not significant"})

    fig = px.scatter(
        plot_df,
        x="best_lfc",
        y="neg_log10_fdr",
        color="significance",
        color_discrete_map={"Significant": "#d62728", "Not significant": "#7f7f7f"},
        hover_name="id",
        hover_data={"best_fdr": ":.2e", "direction": True, "best_lfc": ":.2f"},
        labels={"best_lfc": "Log2 Fold Change (LFC)", "neg_log10_fdr": "-log10(FDR)"},
        title="Volcano Plot",
    )

    fig.add_hline(y=-np.log10(fdr_threshold), line_dash="dash", line_color="gray")
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    fig.update_layout(height=600)

    return fig


def compute_gini_index(values: pd.Series, log_scale: bool = False) -> float:
    """Compute the Gini index of a distribution of read counts.

    0 means perfectly even distribution across sgRNAs; 1 means
    maximally uneven (all reads concentrated in one sgRNA).

    Args:
        values: a Series of non-negative counts (e.g. read counts per sgRNA).
        log_scale: if True, compute the index on log2(count + 1), which is
            how the MAGeCK countsummary reports its GiniIndex. Raw counts
            (the default) are far more sensitive to a few dominant sgRNAs.

    Returns:
        The Gini index as a float between 0 and 1 (NaN if all counts are zero).
    """
    data = values.to_numpy(dtype=float)
    if log_scale:
        data = np.log2(data + 1)

    sorted_values = np.sort(data)
    n = len(sorted_values)
    total = sorted_values.sum()
    if total == 0:
        return float("nan")

    weighted_sum = np.sum(np.arange(1, n + 1) * sorted_values)
    return float((2 * weighted_sum - (n + 1) * total) / (n * total))

def build_read_distribution_plot(df: pd.DataFrame, count_column: str) -> go.Figure:
    """Build a histogram of read counts per sgRNA.

    Args:
        df: DataFrame produced by app.parser.load_sgrna_summary.
        count_column: name of the column holding read counts to plot.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    fig = px.histogram(
        df,
        x=count_column,
        nbins=50,
        labels={count_column: "Read count per sgRNA"},
        title=f"sgRNA Read Distribution ({count_column})",
    )
    fig.update_layout(height=400, showlegend=False)
    return fig


def build_mapping_summary_plot(df: pd.DataFrame) -> go.Figure:
    """Build a grouped bar chart of mapping percentage per sample.

    Args:
        df: DataFrame produced by app.parser.load_countsummary.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    fig = px.bar(
        df,
        x="Label",
        y="Percentage",
        labels={"Label": "Sample", "Percentage": "Mapped reads (%)"},
        title="Mapping Rate per Sample",
        text_auto=".1%",
    )
    fig.update_yaxes(tickformat=".0%", range=[0, 1])
    fig.add_hline(y=0.8, line_dash="dash", line_color="red", annotation_text="80% threshold")
    fig.update_layout(height=400)
    return fig


def build_replicate_correlation_plot(
    df: pd.DataFrame,
    x_column: str,
    y_column: str,
    log_transform: bool = True,
) -> go.Figure:
    """Build a scatter plot comparing read counts between two replicate samples.

    Args:
        df: DataFrame with per-sgRNA count columns (e.g. from a raw
            count table with one column per sample).
        x_column: name of the column for the first replicate.
        y_column: name of the column for the second replicate.
        log_transform: if True (default), apply log2(count + 1) to both
            columns before plotting and computing correlation. This is
            standard practice for count data, which is heavily
            right-skewed and would otherwise be dominated by a few
            highly abundant sgRNAs.

    Returns:
        A Plotly Figure ready to be rendered, including the Pearson
        correlation coefficient in the title.
    """
    plot_df = df.copy()
    x_values, y_values = plot_df[x_column], plot_df[y_column]

    if log_transform:
        x_values = np.log2(x_values + 1)
        y_values = np.log2(y_values + 1)
        plot_df["_x"] = x_values
        plot_df["_y"] = y_values
        x_label, y_label = f"log2({x_column} + 1)", f"log2({y_column} + 1)"
        x_col, y_col = "_x", "_y"
    else:
        x_label, y_label = x_column, y_column
        x_col, y_col = x_column, y_column

    correlation = x_values.corr(y_values)

    fig = px.scatter(
        plot_df,
        x=x_col,
        y=y_col,
        hover_name="sgrna" if "sgrna" in df.columns else None,
        labels={x_col: x_label, y_col: y_label},
        title=f"Replicate Correlation: {x_column} vs {y_column} (Pearson r = {correlation:.3f})",
    )

    # Manual linear trendline (avoids the statsmodels dependency)
    slope, intercept = np.polyfit(x_values, y_values, deg=1)
    x_range = np.array([x_values.min(), x_values.max()])
    y_range = slope * x_range + intercept
    fig.add_trace(
        go.Scatter(
            x=x_range,
            y=y_range,
            mode="lines",
            name="Trendline",
            line=dict(color="gray", dash="dash"),
        )
    )

    fig.update_layout(height=500)
    return fig


def figure_to_png_bytes(fig: go.Figure) -> bytes:
    """Convert a Plotly figure to PNG bytes, for use with st.download_button.

    Args:
        fig: the Plotly figure to export.

    Returns:
        PNG image data as bytes.
    """
    return fig.to_image(format="png", scale=2)


def compare_gene_summaries(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    fdr_threshold: float = 0.05,
    label1: str = "Screen 1",
    label2: str = "Screen 2",
) -> pd.DataFrame:
    """Merge two gene_summary DataFrames for cross-screen comparison.

    Args:
        df1: DataFrame from app.parser.load_gene_summary (first screen).
        df2: DataFrame from app.parser.load_gene_summary (second screen).
        fdr_threshold: FDR cutoff used to flag genes as significant.
        label1: display name for the first screen's columns.
        label2: display name for the second screen's columns.

    Returns:
        A merged DataFrame with one row per gene present in either
        screen, with per-screen LFC/FDR columns and a "concordance"
        column ("Both", "<label1> only", "<label2> only", "Neither").
    """
    merged = df1[["id", "best_lfc", "best_fdr"]].merge(
        df2[["id", "best_lfc", "best_fdr"]],
        on="id",
        how="outer",
        suffixes=(f"_{label1}", f"_{label2}"),
    )

    sig1 = merged[f"best_fdr_{label1}"] < fdr_threshold
    sig2 = merged[f"best_fdr_{label2}"] < fdr_threshold

    conditions = [sig1 & sig2, sig1 & ~sig2, ~sig1 & sig2]
    choices = ["Both", f"{label1} only", f"{label2} only"]
    merged["concordance"] = np.select(conditions, choices, default="Neither")

    return merged


def build_screen_comparison_plot(
    merged_df: pd.DataFrame,
    label1: str = "Screen 1",
    label2: str = "Screen 2",
) -> go.Figure:
    """Build a scatter plot comparing gene LFCs between two screens.

    Args:
        merged_df: output of compare_gene_summaries.
        label1: display name used for the first screen's axis.
        label2: display name used for the second screen's axis.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    fig = px.scatter(
        merged_df,
        x=f"best_lfc_{label1}",
        y=f"best_lfc_{label2}",
        color="concordance",
        hover_name="id",
        labels={f"best_lfc_{label1}": f"LFC ({label1})", f"best_lfc_{label2}": f"LFC ({label2})"},
        title=f"Gene LFC Comparison: {label1} vs {label2}",
        color_discrete_map={
            "Both": "#2ca02c",
            f"{label1} only": "#1f77b4",
            f"{label2} only": "#ff7f0e",
            "Neither": "#cccccc",
        },
    )
    fig.add_shape(
        type="line",
        x0=merged_df[f"best_lfc_{label1}"].min(),
        y0=merged_df[f"best_lfc_{label1}"].min(),
        x1=merged_df[f"best_lfc_{label1}"].max(),
        y1=merged_df[f"best_lfc_{label1}"].max(),
        line=dict(color="gray", dash="dot"),
    )
    fig.update_layout(height=600)
    return fig


def build_zerocounts_plot(df: pd.DataFrame) -> go.Figure:
    """Build a bar chart of sgRNAs with zero reads per sample (log10 scale).

    Args:
        df: DataFrame produced by app.parser.load_countsummary.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    plot_df = df.copy()
    plot_df["log10_zerocounts"] = np.log10(plot_df["Zerocounts"].clip(lower=1))

    fig = px.bar(
        plot_df,
        x="Label",
        y="log10_zerocounts",
        labels={"Label": "Sample", "log10_zerocounts": "log10(sgRNAs with zero reads)"},
        title="sgRNAs with Zero Read Count per Sample",
    )
    fig.update_layout(height=400, showlegend=False)
    return fig


def build_correlation_heatmap(
    df: pd.DataFrame,
    sample_columns: list[str],
    log_transform: bool = True,
) -> go.Figure:
    """Build a clustered-look correlation heatmap across multiple samples.

    Args:
        df: DataFrame with one row per sgRNA and one column per sample.
        sample_columns: names of the numeric columns to correlate.
        log_transform: if True (default), apply log2(count + 1) to each
            column before computing correlations. See
            build_replicate_correlation_plot for rationale.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    values = df[sample_columns]
    if log_transform:
        values = np.log2(values + 1)

    corr_matrix = values.corr()

    fig = px.imshow(
        corr_matrix,
        text_auto=".2f",
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        labels=dict(color="Pearson r"),
        title="Sample Correlation Heatmap" + (" (log2 counts)" if log_transform else " (raw counts)"),
    )
    fig.update_layout(height=500)
    return fig

def build_pca_plot(
    df: pd.DataFrame,
    sample_columns: list[str],
    log_transform: bool = True,
) -> go.Figure:
    """Build a 2D PCA scatter plot of samples based on their sgRNA counts.

    Args:
        df: DataFrame with one row per sgRNA and one column per sample.
        sample_columns: names of the numeric columns to use as features.
        log_transform: if True (default), apply log2(count + 1) before
            PCA. See build_replicate_correlation_plot for rationale;
            without this, PCA is dominated by a handful of very
            highly abundant sgRNAs.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    values = df[sample_columns]
    if log_transform:
        values = np.log2(values + 1)

    data = values.T.to_numpy()

    data_centered = data - data.mean(axis=0)
    stds = data.std(axis=0)
    stds[stds == 0] = 1
    data_scaled = data_centered / stds

    u, s, vt = np.linalg.svd(data_scaled, full_matrices=False)
    scores = u * s
    explained_variance = (s**2) / np.sum(s**2)

    pca_df = pd.DataFrame(
        {
            "sample": sample_columns,
            "PC1": scores[:, 0],
            "PC2": scores[:, 1] if scores.shape[1] > 1 else 0,
        }
    )

    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        text="sample",
        labels={
            "PC1": f"PC1 ({explained_variance[0]:.1%})",
            "PC2": f"PC2 ({explained_variance[1]:.1%})" if len(explained_variance) > 1 else "PC2",
        },
        title="PCA of Samples" + (" (log2 counts)" if log_transform else " (raw counts)"),
    )
    fig.update_traces(textposition="top center", marker=dict(size=12))
    fig.update_layout(height=500)
    return fig

def build_beta_comparison_plot(
    df: pd.DataFrame,
    condition_x: str,
    condition_y: str,
    fdr_threshold: float = 0.05,
    top_n: int = 10,
    fdr_suffix: str = "fdr",
) -> go.Figure:
    """Compare MLE beta scores between two conditions, gene by gene.

    Args:
        df: DataFrame produced by app.parser.load_mle_gene_summary.
        condition_x: condition plotted on the x-axis (e.g. "Control").
        condition_y: condition plotted on the y-axis (e.g. "PLX").
        fdr_threshold: FDR cutoff used to flag genes as significant.
        top_n: number of genes (lowest combined FDR) to label on the plot.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    beta_x, fdr_x = f"{condition_x}|beta", f"{condition_x}|{fdr_suffix}"
    beta_y, fdr_y = f"{condition_y}|beta", f"{condition_y}|{fdr_suffix}"

    plot_df = df.copy()
    sig_x = plot_df[fdr_x] < fdr_threshold
    sig_y = plot_df[fdr_y] < fdr_threshold

    conditions = [sig_x & sig_y, sig_x & ~sig_y, ~sig_x & sig_y]
    choices = ["Both", f"{condition_x} only", f"{condition_y} only"]
    plot_df["significance"] = np.select(conditions, choices, default="Neither")

    fig = px.scatter(
        plot_df,
        x=beta_x,
        y=beta_y,
        color="significance",
        color_discrete_map={
            "Both": "#d62728",
            f"{condition_x} only": "#1f77b4",
            f"{condition_y} only": "#ff7f0e",
            "Neither": "#7f7f7f",
        },
        hover_name="Gene",
        labels={beta_x: f"Beta score ({condition_x})", beta_y: f"Beta score ({condition_y})"},
        title=f"MLE Beta Score Comparison: {condition_x} vs {condition_y}",
    )

    axis_min = min(plot_df[beta_x].min(), plot_df[beta_y].min())
    axis_max = max(plot_df[beta_x].max(), plot_df[beta_y].max())
    fig.add_shape(
        type="line",
        x0=axis_min,
        y0=axis_min,
        x1=axis_max,
        y1=axis_max,
        line=dict(color="gray", dash="dash"),
    )

    plot_df["combined_fdr"] = plot_df[[fdr_x, fdr_y]].min(axis=1)
    top_genes = plot_df.nsmallest(top_n, "combined_fdr")
    fig.add_trace(
        go.Scatter(
            x=top_genes[beta_x],
            y=top_genes[beta_y],
            mode="text",
            text=top_genes["Gene"],
            textposition="top center",
            showlegend=False,
        )
    )

    fig.update_layout(height=600)
    return fig

def build_gc_content_plot(gc_df: pd.DataFrame, sample_name: str = "Sample", nbins: int = 50) -> go.Figure:
    """Build a density plot of per-read GC content.

    Bins the data in Python first (instead of letting Plotly send every
    raw read value to the browser), since GC content is often computed
    over millions of reads for large FASTQ files — sending the raw
    values would exceed Streamlit's message size limit.

    Args:
        gc_df: DataFrame produced by app.fastq_qc.compute_gc_distribution.
        sample_name: label used in the plot title.
        nbins: number of histogram bins.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    counts, bin_edges = np.histogram(gc_df["gc_percent"], bins=nbins, range=(0, 100), density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    fig = px.bar(
        x=bin_centers,
        y=counts,
        labels={"x": "mean % GC content", "y": "density"},
        title=f"GC Content Distribution ({sample_name})",
    )
    fig.update_traces(marker_line_width=0)
    fig.update_layout(height=400, showlegend=False, bargap=0)
    return fig

def build_quality_profile_plot(quality_df: pd.DataFrame, sample_name: str = "Sample") -> go.Figure:
    """Build a line plot of median base quality by read position.

    Args:
        quality_df: DataFrame produced by app.fastq_qc.compute_quality_profile.
        sample_name: label used in the plot title.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    fig = px.line(
        quality_df,
        x="position",
        y="median_quality",
        labels={"position": "base position", "median_quality": "median quality"},
        title=f"Base Quality Distribution ({sample_name})",
    )
    fig.update_layout(height=400, showlegend=False, yaxis_range=[0, 45])
    return fig

def build_beta_score_plot(
    df: pd.DataFrame,
    condition: str,
    fdr_threshold: float = 0.05,
    fdr_suffix: str = "fdr",
) -> go.Figure:
    """Build a volcano-style plot of MAGeCK MLE beta scores for one condition.

    Args:
        df: DataFrame produced by app.parser.load_mle_gene_summary.
        condition: name of the condition to plot (e.g. "tumor").
        fdr_threshold: FDR cutoff used to color genes as significant.
        fdr_suffix: FDR column to use: "fdr" (permutation) or "wald-fdr".

    Returns:
        A Plotly Figure ready to be rendered.
    """
    beta_col = f"{condition}|beta"
    fdr_col = f"{condition}|{fdr_suffix}"

    plot_df = df.copy()
    fdr_values = plot_df[fdr_col]

    # An FDR of exactly 0 would map to an extreme y value and squash every other gene,
    # so floor it one order of magnitude below the smallest non-zero FDR in the file.
    smallest_nonzero = fdr_values[fdr_values > 0].min()
    floor = smallest_nonzero / 10 if pd.notna(smallest_nonzero) else 1e-300
    plot_df["neg_log10_fdr"] = -np.log10(fdr_values.clip(lower=floor))

    is_significant = fdr_values < fdr_threshold
    plot_df["significance"] = is_significant.map({True: "Significant", False: "Not significant"})

    fig = px.scatter(
        plot_df,
        x=beta_col,
        y="neg_log10_fdr",
        color="significance",
        color_discrete_map={"Significant": "#d62728", "Not significant": "#7f7f7f"},
        hover_name="Gene",
        hover_data={fdr_col: ":.2e"},
        labels={beta_col: f"Beta score ({condition})", "neg_log10_fdr": "-log10(FDR)"},
        title=f"MLE Beta Score Plot: {condition}",
    )
    fig.add_hline(y=-np.log10(fdr_threshold), line_dash="dash", line_color="gray")
    fig.update_layout(height=600)
    return fig