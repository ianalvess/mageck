import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from app.library_check import (
    filter_observed,
    observed_mask,
    summarize_observed,
    write_observed_count_table,
)

counts = pd.DataFrame(
    {
        "sgRNA": ["HGLibA_1", "HGLibA_2", "HGLibA_3", "HGLibB_1", "HGLibB_2"],
        "Gene": ["G1", "G1", "G2", "G1", "G2"],
        "control1": [10, 0, 5, 0, 0],
        "tumor1": [0, 0, 7, 0, 0],
    }
)
samples = ["control1", "tumor1"]

# HGLibA_2 and both HGLibB sgRNAs have no reads in any sample
assert observed_mask(counts, samples).tolist() == [True, False, True, False, False]

summary = summarize_observed(counts, samples)
print(summary["by_prefix"])
assert summary["total"] == 5
assert summary["unobserved"] == 3
assert summary["by_prefix"].loc["HGLibA", "unobserved"] == 1
assert summary["by_prefix"].loc["HGLibB", "observed"] == 0
assert summary["by_prefix"].loc["HGLibB", "unobserved"] == 2

filtered = filter_observed(counts, samples)
assert filtered["sgRNA"].tolist() == ["HGLibA_1", "HGLibA_3"]
print("OK: unobserved sgRNAs are detected, summarized by prefix and filtered")

with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "counts.txt"
    target = Path(tmp) / "out" / "counts.observed.txt"
    counts.to_csv(source, sep="\t", index=False)

    result = write_observed_count_table(source, target)
    written = pd.read_csv(target, sep="\t")

    assert result == {"kept": 2, "removed": 3}
    assert written["sgRNA"].tolist() == ["HGLibA_1", "HGLibA_3"]
    assert written.columns.tolist() == ["sgRNA", "Gene", "control1", "tumor1"]
print("OK: filtered count table written with the expected sgRNAs and columns")

from app.library_check import count_sample_columns

count_table = pd.DataFrame(
    {"sgRNA": ["a", "b"], "Gene": ["G1", "G2"], "control1": [1, 2], "tumor1": [3, 4]}
)
assert count_sample_columns(count_table) == ["control1", "tumor1"]

gene_summary = pd.DataFrame({"id": ["G1"], "num": [3], "neg|score": [0.1], "neg|fdr": [0.5]})
assert count_sample_columns(gene_summary) == []

sgrna_summary = pd.DataFrame({"sgrna": ["a"], "Gene": ["G1"], "LFC": [1.0], "FDR": [0.1]})
assert count_sample_columns(sgrna_summary) == []

count_summary = pd.DataFrame({"Label": ["c1"], "Reads": [10], "GiniIndex": [0.1], "Zerocounts": [5]})
assert count_sample_columns(count_summary) == []

assert count_sample_columns(count_table.assign(tumor1=[-1, 4])) == ["control1"]
print("OK: only per-sample count columns are accepted")