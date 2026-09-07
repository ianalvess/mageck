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

    # Avoid log10(0) errors when FDR is exactly 0
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

    # Manual linear trendline (avoids the statsmodels dependency)
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
        # Only embed the Plotly JS library once, in the first chart
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