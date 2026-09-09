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
        labels={"best_lfc": "Log Fold Change", "neg_log10_fdr": "-log10(FDR)"},
        title="Volcano Plot",
    )

    fig.add_hline(y=-np.log10(fdr_threshold), line_dash="dash", line_color="gray")
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    fig.update_layout(height=600)

    return fig


def compute_gini_index(values: pd.Series) -> float:
    """Compute the Gini index of a distribution of read counts.

    0 means perfectly even distribution across sgRNAs; 1 means
    maximally uneven (all reads concentrated in one sgRNA).

    Args:
        values: a Series of non-negative counts (e.g. read counts per sgRNA).

    Returns:
        The Gini index as a float between 0 and 1.
    """
    sorted_values = np.sort(values.to_numpy())
    n = len(sorted_values)
    cumulative = np.cumsum(sorted_values)
    return (2 * np.sum((np.arange(1, n + 1)) * sorted_values) - (n + 1) * cumulative[-1]) / (
        n * cumulative[-1]
    )


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
) -> go.Figure:
    """Build a scatter plot comparing read counts between two replicate samples.

    Args:
        df: DataFrame with per-sgRNA count columns (e.g. from a raw
            count table with one column per sample).
        x_column: name of the column for the first replicate.
        y_column: name of the column for the second replicate.

    Returns:
        A Plotly Figure ready to be rendered, including the Pearson
        correlation coefficient in the title.
    """
    correlation = df[x_column].corr(df[y_column])

    fig = px.scatter(
        df,
        x=x_column,
        y=y_column,
        hover_name="sgrna" if "sgrna" in df.columns else None,
        labels={x_column: x_column, y_column: y_column},
        title=f"Replicate Correlation: {x_column} vs {y_column} (r = {correlation:.3f})",
    )

    slope, intercept = np.polyfit(df[x_column], df[y_column], deg=1)
    x_range = np.array([df[x_column].min(), df[x_column].max()])
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


def build_html_report(figures: dict[str, go.Figure], title: str = "MAGeCK Dashboard Report") -> str:
    """Combine multiple Plotly figures into a single standalone HTML report.

    Args:
        figures: mapping of section title (e.g. "Volcano Plot") to the
            corresponding Plotly Figure. Only non-None entries are included.
        title: title shown at the top of the report.

    Returns:
        A complete HTML document as a string, with Plotly's JS embedded
        so it can be opened offline in any browser.
    """
    sections = []
    for i, (section_title, fig) in enumerate(figures.items()):
        if fig is None:
            continue
        include_js = "cdn" if i == 0 else False
        chart_html = fig.to_html(full_html=False, include_plotlyjs=include_js)
        sections.append(f"<h2>{section_title}</h2>\n{chart_html}")

    body = "\n<hr>\n".join(sections)

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
    body {{ font-family: sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; }}
    h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
    h2 {{ margin-top: 40px; color: #444; }}
</style>
</head>
<body>
<h1>{title}</h1>
{body}
</body>
</html>"""


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


def build_correlation_heatmap(df: pd.DataFrame, sample_columns: list[str]) -> go.Figure:
    """Build a clustered-look correlation heatmap across multiple samples.

    Args:
        df: DataFrame with one row per sgRNA and one column per sample.
        sample_columns: names of the numeric columns to correlate.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    corr_matrix = df[sample_columns].corr()

    fig = px.imshow(
        corr_matrix,
        text_auto=".2f",
        color_continuous_scale="RdBu_r",
        zmin=-1,
        zmax=1,
        labels=dict(color="Pearson r"),
        title="Sample Correlation Heatmap",
    )
    fig.update_layout(height=500)
    return fig


def build_pca_plot(df: pd.DataFrame, sample_columns: list[str]) -> go.Figure:
    """Build a 2D PCA scatter plot of samples based on their sgRNA counts.

    Args:
        df: DataFrame with one row per sgRNA and one column per sample.
        sample_columns: names of the numeric columns to use as features.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    data = df[sample_columns].T.to_numpy()

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
        title="PCA of Samples",
    )
    fig.update_traces(textposition="top center", marker=dict(size=12))
    fig.update_layout(height=500)
    return fig


def build_beta_score_plot(
    df: pd.DataFrame,
    condition: str,
    fdr_threshold: float = 0.05,
) -> go.Figure:
    """Build a volcano-style plot of MAGeCK MLE beta scores for one condition.

    Args:
        df: DataFrame produced by app.parser.load_mle_gene_summary.
        condition: name of the condition to plot (e.g. "PLX").
        fdr_threshold: FDR cutoff used to color genes as significant.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    beta_col = f"{condition}|beta"
    fdr_col = f"{condition}|fdr"

    plot_df = df.copy()
    plot_df["neg_log10_fdr"] = -np.log10(plot_df[fdr_col].clip(lower=1e-300))
    is_significant = plot_df[fdr_col] < fdr_threshold
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


def build_gc_content_plot(gc_df: pd.DataFrame, sample_name: str = "Sample") -> go.Figure:
    """Build a density plot of per-read GC content.

    Args:
        gc_df: DataFrame produced by app.fastq_qc.compute_gc_distribution.
        sample_name: label used in the plot title.

    Returns:
        A Plotly Figure ready to be rendered.
    """
    fig = px.histogram(
        gc_df,
        x="gc_percent",
        histnorm="probability density",
        nbins=50,
        labels={"gc_percent": "mean % GC content"},
        title=f"GC Content Distribution ({sample_name})",
    )
    fig.update_traces(marker_line_width=0)
    fig.update_layout(height=400, showlegend=False)
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