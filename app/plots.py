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