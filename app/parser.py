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
    - direction: "neg" or "pos", chosen by comparing the two FDRs of the gene, not by the
      sign of any LFC. The smaller FDR wins and exact ties go to "neg". A missing FDR
      never wins, and the value is NaN when neither FDR is available. Every gene with at
      least one FDR gets a direction, including non-significant genes, so for non-hits
      it is not evidence of a biological effect.
    - best_fdr: the FDR of the chosen direction, copied unchanged. MAGeCK calculates the
      neg (depletion) and pos (enrichment) FDRs separately as one-sided tests, so
      best_fdr is the smaller of two separately calculated one-sided FDRs. It is not a
      separately corrected two-sided FDR, and it does not account for having picked
      the better of the two directions.
    - best_lfc: the LFC of the same direction (neg|lfc or pos|lfc), so its sign can
      disagree with the direction.

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

    # Pick whichever direction (depletion vs enrichment) is more significant per gene.
    # A comparison with a missing (NaN) FDR is always False, so a gene with a valid neg FDR
    # and a missing pos FDR must be assigned to neg explicitly; a gene with no FDR at all
    # has no direction.
    neg_fdr, pos_fdr = df["neg|fdr"], df["pos|fdr"]
    neg_is_better = (neg_fdr <= pos_fdr) | (neg_fdr.notna() & pos_fdr.isna())

    df["direction"] = neg_is_better.map({True: "neg", False: "pos"}).where(neg_fdr.notna() | pos_fdr.notna())
    df["best_fdr"] = neg_fdr.where(neg_is_better, pos_fdr)
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


def load_countsummary(path: str | Path) -> pd.DataFrame:
    """Load and clean a MAGeCK countsummary.txt file.

    This file summarizes read mapping quality per sample, produced
    by `mageck count` (runs before `mageck test`).

    Args:
        path: path to the countsummary.txt file.

    Returns:
        Cleaned DataFrame with the original columns.
    """
    df = pd.read_csv(path, sep="\t", na_values=["NA", "na", ""])

    required = {"Label", "Reads", "Mapped", "Percentage", "Zerocounts", "GiniIndex"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"countsummary file is missing expected columns: {missing}")

    return df


def load_mle_gene_summary(path: str | Path) -> tuple[pd.DataFrame, list[str]]:
    """Load and clean a MAGeCK MLE gene_summary.txt file.

    Unlike `mageck test`, `mageck mle` output can contain one or more
    conditions, each with its own <condition>|beta, <condition>|z,
    <condition>|p-value, <condition>|fdr columns.

    Args:
        path: path to the gene_summary.txt file produced by `mageck mle`.

    Returns:
        A tuple of (DataFrame, list of condition names found in the file).
    """
    df = pd.read_csv(path, sep="\t", na_values=["NA", "na", ""])

    if "Gene" not in df.columns:
        raise ValueError("MLE gene_summary file is missing the expected 'Gene' column")

    conditions = sorted({col.split("|")[0] for col in df.columns if col.endswith("|beta")})

    if not conditions:
        raise ValueError("No <condition>|beta columns found in this MLE gene_summary file")

    return df, conditions