"""
Helpers to detect sgRNAs that were never observed in any sample.

A library file often contains more sgRNAs than were actually used in a
screen (for example, one half of a split library). Those sgRNAs have zero
reads in every sample, which inflates zero-count statistics and distorts
normalization, correlation and PCA.
"""

from pathlib import Path

import pandas as pd


def observed_mask(counts_df: pd.DataFrame, sample_columns: list[str]) -> pd.Series:
    """Return a boolean mask: True for sgRNAs with at least one read in any sample."""
    return counts_df[sample_columns].sum(axis=1) > 0


def summarize_observed(
    counts_df: pd.DataFrame,
    sample_columns: list[str],
    id_column: str = "sgRNA",
) -> dict:
    """Summarize how many sgRNAs were never observed, grouped by ID prefix.

    The prefix is the part of the sgRNA ID before the first underscore
    (e.g. "HGLibA" in "HGLibA_00001"), which usually identifies the
    library half or sub-pool.

    Returns:
        A dict with total, unobserved, unobserved_fraction and a by_prefix
        DataFrame (columns: total, observed, unobserved).
    """
    mask = observed_mask(counts_df, sample_columns)
    prefixes = counts_df[id_column].astype(str).str.split("_").str[0]

    by_prefix = (
        pd.DataFrame({"prefix": prefixes, "observed": mask})
        .groupby("prefix")["observed"]
        .agg(total="size", observed="sum")
    )
    by_prefix["unobserved"] = by_prefix["total"] - by_prefix["observed"]

    return {
        "total": int(len(counts_df)),
        "unobserved": int((~mask).sum()),
        "unobserved_fraction": float((~mask).mean()),
        "by_prefix": by_prefix,
    }


def filter_observed(counts_df: pd.DataFrame, sample_columns: list[str]) -> pd.DataFrame:
    """Return only the sgRNAs observed in at least one sample."""
    return counts_df[observed_mask(counts_df, sample_columns)].copy()


def write_observed_count_table(count_path: str | Path, output_path: str | Path) -> dict:
    """Write a copy of a MAGeCK count table keeping only observed sgRNAs.

    Args:
        count_path: path to a tab-separated count table (sgRNA, Gene, samples...).
        output_path: where to write the filtered table.

    Returns:
        A dict with the number of sgRNAs kept and removed.

    Raises:
        ValueError: if no per-sample count columns are found in the table.
    """
    # Read the identifier columns (sgRNA, Gene) as text so numeric-looking IDs are
    # neither mistaken for samples nor altered (e.g. leading zeros) when written back.
    header = pd.read_csv(count_path, sep="\t", nrows=0).columns.tolist()
    counts_df = pd.read_csv(count_path, sep="\t", dtype={column: str for column in header[:2]})
    sample_columns = count_sample_columns(counts_df)
    if not sample_columns:
        raise ValueError("No per-sample count columns found in this count table")
    filtered = filter_observed(counts_df, sample_columns)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(output_path, sep="\t", index=False)

    return {"kept": len(filtered), "removed": len(counts_df) - len(filtered)}

SUMMARY_FILE_MARKERS = {
    "lfc", "fdr", "score", "p.low", "p.high", "p.twosided",
    "control_count", "treatment_count", "control_mean", "treat_mean",
    "giniindex", "zerocounts", "negselqc",
}


def count_sample_columns(df: pd.DataFrame) -> list[str]:
    """Return the per-sample count columns of a MAGeCK count table.

    Tables that look like a MAGeCK gene/sgRNA summary or a countsummary (columns
    such as FDR, LFC, score, GiniIndex or "neg|fdr") are rejected, so statistics
    are never offered as samples. The first two columns (sgRNA and Gene) are
    identifiers and are never returned, even when their values look numeric.

    Returns:
        The numeric, non-negative sample column names, or an empty list if the
        table is not a count table.
    """
    names = [str(column).lower() for column in df.columns]
    if any("|" in name or name in SUMMARY_FILE_MARKERS for name in names):
        return []

    numeric = df.iloc[:, 2:].select_dtypes(include="number")
    return [column for column in numeric.columns if (numeric[column].dropna() >= 0).all()]