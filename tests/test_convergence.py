import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from app.convergence import find_convergent_sgrnas, summarize_gene_convergence

counts = pd.DataFrame(
    {
        "sgRNA": [f"s{i}" for i in range(1, 13)],
        "Gene": ["G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4", "G5", "G5", "G6", "G6"],
        "control1": [100, 100, 100, 100, 50, 100, 100, 100, 100, 100, 100, 100],
        "control2": [100, 100, 100, 100, 50, 100, 100, 100, 100, 100, 100, 100],
        "tumor1": [5000, 4000, 8000, 90, 2000, 100, 100, 100, 150, 130, 100, 100],
        "tumor2": [6000, 3000, 20, 100, 2500, 120, 100, 100, 130, 150, 100, 100],
    }
)
controls = ["control1", "control2"]
replicates = ["tumor1", "tumor2"]

# s3 expands only in tumor1 (a one-replicate jackpot) and must NOT qualify
result = find_convergent_sgrnas(counts, controls, replicates, normalize_depth=False)
print(result[["sgRNA", "Gene", "fold_tumor1", "fold_tumor2"]])
assert result["sgRNA"].tolist() == ["s1", "s2", "s5"]
print("OK: only sgRNAs enriched in every replicate qualify")

summary = summarize_gene_convergence(result, counts)
assert summary["n_sgrnas"] == 3
assert summary["observed_pairs"] == 1  # s1 and s2 both belong to G1
assert abs(summary["expected_pairs"] - 3 / 11) < 1e-9
assert summary["genes"]["Gene"].tolist() == ["G1"]
assert summary["genes"]["n_convergent"].tolist() == [2]
assert summary["genes"]["n_in_library"].tolist() == [2]
print("OK: gene-level convergence summarized against the chance expectation")

empty_column = counts.assign(tumor2=0)
try:
    find_convergent_sgrnas(empty_column, controls, replicates)
    raise AssertionError("Expected ValueError for a column with no reads")
except ValueError:
    print("OK: a column with no reads raises ValueError")