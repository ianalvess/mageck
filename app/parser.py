"""
Parsers for MAGeCK output files.

MAGeCK produces two main tab-separated files after `mageck test`:
- <prefix>.gene_summary.txt : one row per gene, with separate
  depletion ("neg") and enrichment ("pos") statistics.
- <prefix>.sgrna_summary.txt : one row per individual sgRNA.

These functions only depend on pandas so they can be tested and
reused outside of the Streamlit app.
"""

from pathlib import Path

import pandas as pd


def load_gene_summary(path: str | Path) -> pd.DataFrame:
    """Load and clean a MAGeCK gene_summary.txt file.

    Adds convenience columns:
    - direction: "neg" or "pos", whichever has the lower FDR for that gene
    - best_fdr / best_lfc: the FDR/LFC value for the chosen direction

    Args:
        path: path to the gene_summary.txt file.

    Returns:
        Cleaned DataFrame with the original columns plus direction,
        best_fdr and best_lfc.
    """
    df = pd.read_csv(path, sep="\t", na_values=["NA", "na", ""])

    required = {"id", "neg|score", "neg|lfc", "neg|fdr", "pos|score", "pos|lfc", "pos|fdr"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"gene_summary file is missing expected columns: {missing}")

    # Pick whichever direction (depletion vs enrichment) is more significant per gene
    neg_is_better = df["neg|fdr"] <= df["pos|fdr"]

    df["direction"] = neg_is_better.map({True: "neg", False: "pos"})
    df["best_fdr"] = df["neg|fdr"].where(neg_is_better, df["pos|fdr"])
    df["best_lfc"] = df["neg|lfc"].where(neg_is_better, df["pos|lfc"])

    return df


def load_sgrna_summary(path: str | Path) -> pd.DataFrame:
    """Load and clean a MAGeCK sgrna_summary.txt file.

    Args:
        path: path to the sgrna_summary.txt file.

    Returns:
        Cleaned DataFrame with numeric columns coerced properly.
    """
    df = pd.read_csv(path, sep="\t", na_values=["NA", "na", ""])

    required = {"sgrna", "Gene", "LFC", "FDR"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"sgrna_summary file is missing expected columns: {missing}")

    return df