"""
Replicate convergence analysis for pooled CRISPR screens.

Gene-level tests can be dominated by a single replicate in which one
sgRNA expanded hugely (a "jackpot" clone). A more robust signal is an
sgRNA enriched in every replicate of a condition, and a gene whose
independent sgRNAs converge on the same result.
"""

from math import comb

import pandas as pd


def find_convergent_sgrnas(
    counts_df: pd.DataFrame,
    control_columns: list[str],
    replicate_columns: list[str],
    min_reads: int = 100,
    min_fold: float = 10.0,
    normalize_depth: bool = True,
    id_column: str = "sgRNA",
    gene_column: str = "Gene",
) -> pd.DataFrame:
    """Find sgRNAs enriched in every replicate relative to the controls.

    An sgRNA qualifies if, in each replicate column, it has at least
    min_reads raw reads and its count is more than min_fold times the mean
    of the control columns (with a pseudocount of 1).

    Args:
        counts_df: count table (one row per sgRNA, one column per sample).
        control_columns: control sample columns.
        replicate_columns: replicate columns of the condition to test.
        min_reads: minimum raw reads required in every replicate.
        min_fold: minimum fold over the control mean required in every replicate.
        normalize_depth: scale each column to the mean library size before
            computing folds.
        id_column: column with sgRNA IDs.
        gene_column: column with gene names.

    Returns:
        The qualifying sgRNAs with their counts and per-replicate folds.
    """
    all_columns = control_columns + replicate_columns
    scaled = counts_df[all_columns].astype(float)

    if normalize_depth:
        totals = scaled.sum()
        if (totals == 0).any():
            raise ValueError("A selected column has no reads, so it cannot be depth-normalized")
        scaled = scaled / (totals / totals.mean())

    control_mean = scaled[control_columns].mean(axis=1)

    qualifies = pd.Series(True, index=counts_df.index)
    folds = {}
    for column in replicate_columns:
        fold = (scaled[column] + 1) / (control_mean + 1)
        folds[f"fold_{column}"] = fold
        qualifies &= (counts_df[column] >= min_reads) & (fold > min_fold)

    result = counts_df.loc[qualifies, [id_column, gene_column] + all_columns].copy()
    for name, fold in folds.items():
        result[name] = fold[qualifies]

    return result.sort_values(gene_column).reset_index(drop=True)


def summarize_gene_convergence(
    convergent_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    gene_column: str = "Gene",
) -> dict:
    """Summarize genes with two or more convergent sgRNAs.

    Compares the observed number of same-gene sgRNA pairs with the number
    expected if the convergent sgRNAs were picked at random from the library.

    Returns:
        A dict with n_sgrnas, observed_pairs, expected_pairs and genes
        (a DataFrame of genes with 2+ convergent sgRNAs).
    """
    n_total = len(counts_df)
    sgrnas_per_gene = n_total / counts_df[gene_column].nunique()

    k = len(convergent_df)
    per_gene = convergent_df[gene_column].value_counts()
    observed_pairs = int(sum(comb(int(n), 2) for n in per_gene))
    expected_pairs = comb(k, 2) * (sgrnas_per_gene - 1) / (n_total - 1) if n_total > 1 else 0.0

    in_library = counts_df[gene_column].value_counts()
    multi = per_gene[per_gene >= 2]
    genes = pd.DataFrame(
        {
            gene_column: multi.index,
            "n_convergent": multi.to_numpy(),
            "n_in_library": [int(in_library[gene]) for gene in multi.index],
        }
    )

    return {
        "n_sgrnas": k,
        "observed_pairs": observed_pairs,
        "expected_pairs": float(expected_pairs),
        "genes": genes,
    }