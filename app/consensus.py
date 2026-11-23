"""
Gene-level consensus across analysis methods.

Combines the RRA results (`mageck test`), the MLE results (`mageck mle`) and the
replicate convergence analysis into one table, with the number of methods that
support each gene at the chosen thresholds. The consensus ranks candidates; it does
not validate them, since all methods share the same replicates.
"""

from functools import reduce

import numpy as np
import pandas as pd


def find_duplicate_genes(
    rra_df: pd.DataFrame | None = None,
    mle_df: pd.DataFrame | None = None,
) -> dict[str, list[str]]:
    """Gene IDs that appear more than once in the RRA or MLE input.

    Returns {"rra": [...], "mle": [...]} with only the sources that have duplicates.
    The convergence table is not checked: it has one row per sgRNA, so a gene is
    expected to repeat there.
    """
    found = {}
    for source, frame, column in (("rra", rra_df, "id"), ("mle", mle_df, "Gene")):
        if frame is None:
            continue
        counts = frame[column].dropna().astype(str).value_counts()
        repeated = sorted(counts[counts > 1].index)
        if repeated:
            found[source] = repeated
    return found


def describe_duplicate_genes(duplicates: dict[str, list[str]], max_listed: int = 5) -> str:
    """One-line description of find_duplicate_genes output, for errors."""
    labels = {"rra": "RRA (mageck test)", "mle": "MLE (mageck mle)"}
    parts = []
    for source, genes in duplicates.items():
        listed = ", ".join(genes[:max_listed])
        extra = f" and {len(genes) - max_listed} more" if len(genes) > max_listed else ""
        parts.append(f"{labels[source]}: {listed}{extra}")
    return "Duplicate gene IDs found. " + "; ".join(parts) + ". Each gene must appear once per input."


def _gene_set(values: pd.Series) -> set[str]:
    return set(values.dropna().astype(str))


def _normalized_gene_set(genes: set[str]) -> set[str]:
    return {gene.strip().casefold() for gene in genes}


def gene_overlap_report(
    rra_df: pd.DataFrame | None = None,
    mle_df: pd.DataFrame | None = None,
    convergent_df: pd.DataFrame | None = None,
    gene_column: str = "Gene",
    min_overlap: float = 0.5,
) -> list[str]:
    """Warnings for inputs whose gene names barely overlap.

    Compares each pair of provided inputs (RRA, MLE, convergent genes) by exact gene
    name. The overlap is shared genes divided by the size of the smaller set; a pair
    below min_overlap gets one warning. Nothing is merged or renamed: when ignoring
    case and surrounding spaces would raise the overlap, the warning says so.
    """
    sets = []
    if rra_df is not None:
        sets.append(("RRA", _gene_set(rra_df["id"])))
    if mle_df is not None:
        sets.append(("MLE", _gene_set(mle_df["Gene"])))
    if convergent_df is not None:
        sets.append(("convergence", _gene_set(convergent_df[gene_column])))

    warnings = []
    for i, (name_a, genes_a) in enumerate(sets):
        for name_b, genes_b in sets[i + 1:]:
            if not genes_a or not genes_b:
                continue
            smaller = min(len(genes_a), len(genes_b))
            shared = len(genes_a & genes_b)
            if shared / smaller >= min_overlap:
                continue
            message = (
                f"Only {shared} of the {smaller} genes in the smaller set ({shared / smaller:.0%}) are shared "
                f"between {name_a} and {name_b}. Gene names are matched exactly, so these inputs may come "
                "from different libraries or use different naming."
            )
            loose = len(_normalized_gene_set(genes_a) & _normalized_gene_set(genes_b))
            if loose > shared:
                message += f" Ignoring case and surrounding spaces would match {loose}."
            warnings.append(message)
    return warnings


def build_consensus_table(
    rra_df: pd.DataFrame | None = None,
    mle_df: pd.DataFrame | None = None,
    convergent_df: pd.DataFrame | None = None,
    mle_condition: str | None = None,
    mle_fdr_suffix: str = "fdr",
    rra_fdr_threshold: float = 0.05,
    mle_fdr_threshold: float = 0.05,
    min_convergent_sgrnas: int = 2,
    gene_column: str = "Gene",
    mle_max_abs_beta: float | None = 10.0,
) -> pd.DataFrame:
    """Cross-reference the analysis methods gene by gene.

    Args:
        rra_df: DataFrame from app.parser.load_gene_summary (mageck test).
        mle_df: DataFrame from app.parser.load_mle_gene_summary (mageck mle).
        convergent_df: sgRNA-level table from app.convergence.find_convergent_sgrnas.
        mle_condition: MLE condition to use (required when mle_df is given).
        mle_fdr_suffix: MLE FDR column to use: "fdr" (permutation) or "wald-fdr".
        rra_fdr_threshold: a gene is supported by RRA if its best FDR is below this.
        mle_fdr_threshold: a gene is supported by MLE if its FDR is below this.
        min_convergent_sgrnas: a gene is supported by convergence if at least this
            many of its sgRNAs are enriched in every replicate.
        gene_column: gene column of convergent_df.
        mle_max_abs_beta: MLE betas above this absolute value are flagged as
            unstable (typical of sgRNAs with near-zero counts in one sample) and
            do not count as support. Use None to disable the cap.

    Returns:
        One row per gene, sorted by number of supporting methods and then by the
        smallest FDR, with the columns of each method that was provided,
        methods_supporting and consensus_direction ("enriched", "depleted",
        "mixed" when the supporting methods disagree, or "" with no support).
    """
    if rra_df is None and mle_df is None and convergent_df is None:
        raise ValueError("Provide at least one of rra_df, mle_df or convergent_df")
    if mle_df is not None and not mle_condition:
        raise ValueError("mle_condition is required when mle_df is provided")
    duplicates = find_duplicate_genes(rra_df, mle_df)
    if duplicates:
        raise ValueError(describe_duplicate_genes(duplicates))

    frames = []
    if rra_df is not None:
        frames.append(
            pd.DataFrame(
                {
                    "Gene": rra_df["id"],
                    "rra_fdr": rra_df["best_fdr"],
                    "rra_lfc": rra_df["best_lfc"],
                    "rra_direction": rra_df["direction"].map({"neg": "depleted", "pos": "enriched"}),
                }
            )
        )
    if mle_df is not None:
        frames.append(
            pd.DataFrame(
                {
                    "Gene": mle_df["Gene"],
                    "mle_beta": mle_df[f"{mle_condition}|beta"],
                    "mle_fdr": mle_df[f"{mle_condition}|{mle_fdr_suffix}"],
                }
            )
        )
    if convergent_df is not None:
        frames.append(
            convergent_df[gene_column].value_counts().rename_axis("Gene").reset_index(name="n_convergent")
        )

    table = reduce(lambda left, right: left.merge(right, on="Gene", how="outer"), frames)

    support = {}
    directions = []
    if rra_df is not None:
        support["rra_support"] = table["rra_fdr"] < rra_fdr_threshold
        directions.append(table["rra_direction"].where(support["rra_support"]))
    if mle_df is not None:
        if mle_max_abs_beta is None:
            table["mle_unstable"] = False
        else:
            table["mle_unstable"] = table["mle_beta"].abs() > mle_max_abs_beta
        support["mle_support"] = (table["mle_fdr"] < mle_fdr_threshold) & ~table["mle_unstable"]
        mle_direction = pd.Series(np.where(table["mle_beta"] > 0, "enriched", "depleted"), index=table.index)
        directions.append(mle_direction.where(support["mle_support"]))
    if convergent_df is not None:
        table["n_convergent"] = table["n_convergent"].fillna(0).astype(int)
        support["convergence_support"] = table["n_convergent"] >= min_convergent_sgrnas
        # Convergence only looks for enrichment, so it can only support "enriched"
        directions.append(pd.Series("enriched", index=table.index).where(support["convergence_support"]))

    for name, values in support.items():
        table[name] = values
    table["methods_supporting"] = table[list(support)].sum(axis=1).astype(int)

    def combine_directions(row) -> str:
        found = {value for value in row if isinstance(value, str)}
        if not found:
            return ""
        return found.pop() if len(found) == 1 else "mixed"

    stacked = pd.concat(directions, axis=1, ignore_index=True)
    table["consensus_direction"] = stacked.apply(combine_directions, axis=1) if len(stacked) else ""

    fdr_columns = [column for column in ("rra_fdr", "mle_fdr") if column in table.columns]
    table["min_fdr"] = table[fdr_columns].min(axis=1) if fdr_columns else np.nan
    table = table.sort_values(["methods_supporting", "min_fdr"], ascending=[False, True])

    ordered = [
        "Gene", "methods_supporting", "consensus_direction",
        "rra_support", "rra_fdr", "rra_lfc", "rra_direction",
        "mle_support", "mle_beta", "mle_fdr", "mle_unstable",
        "convergence_support", "n_convergent", "min_fdr",
    ]
    return table[[column for column in ordered if column in table.columns]].reset_index(drop=True)