import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from app.consensus import build_consensus_table

rra = pd.DataFrame(
    {
        "id": ["G1", "G2", "G3", "G4"],
        "best_fdr": [0.01, 0.2, 0.01, 0.5],
        "best_lfc": [3.0, 1.0, -2.0, 0.1],
        "direction": ["pos", "pos", "neg", "pos"],
    }
)
mle = pd.DataFrame(
    {
        "Gene": ["G1", "G2", "G3", "G5"],
        "T|beta": [5.0, 4.0, 1.0, -2.0],
        "T|fdr": [0.01, 0.01, 0.02, 0.03],
        "T|wald-fdr": [0.9, 0.9, 0.9, 0.9],
    }
)
convergent = pd.DataFrame({"sgRNA": ["a", "b", "c"], "Gene": ["G1", "G1", "G2"]})

table = build_consensus_table(rra_df=rra, mle_df=mle, convergent_df=convergent, mle_condition="T")
print(table[["Gene", "methods_supporting", "consensus_direction", "n_convergent"]])

assert table["Gene"].tolist() == ["G1", "G3", "G2", "G5", "G4"]
rows = table.set_index("Gene")
assert rows.loc["G1", "methods_supporting"] == 3
assert rows.loc["G1", "consensus_direction"] == "enriched"
assert rows.loc["G3", "consensus_direction"] == "mixed"  # RRA says depleted, MLE says enriched
assert rows.loc["G5", "consensus_direction"] == "depleted"  # only MLE supports it
assert rows.loc["G4", "methods_supporting"] == 0 and rows.loc["G4", "consensus_direction"] == ""
assert rows.loc["G2", "n_convergent"] == 1 and rows.loc["G5", "n_convergent"] == 0
print("OK: consensus counts supporting methods and flags direction conflicts")

wald_table = build_consensus_table(mle_df=mle, mle_condition="T", mle_fdr_suffix="wald-fdr")
assert wald_table["methods_supporting"].sum() == 0
print("OK: the MLE FDR column follows mle_fdr_suffix")

for bad_call in (lambda: build_consensus_table(), lambda: build_consensus_table(mle_df=mle)):
    try:
        bad_call()
        raise AssertionError("Expected ValueError")
    except ValueError:
        pass
print("OK: missing inputs raise ValueError")

unstable = pd.DataFrame({"Gene": ["U1", "U2"], "T|beta": [20.0, 2.0], "T|fdr": [0.001, 0.001]})
capped = build_consensus_table(mle_df=unstable, mle_condition="T").set_index("Gene")
assert bool(capped.loc["U1", "mle_unstable"]) and capped.loc["U1", "methods_supporting"] == 0
assert not bool(capped.loc["U2", "mle_unstable"]) and capped.loc["U2", "methods_supporting"] == 1
uncapped = build_consensus_table(mle_df=unstable, mle_condition="T", mle_max_abs_beta=None).set_index("Gene")
assert uncapped.loc["U1", "methods_supporting"] == 1
print("OK: unstable MLE betas are flagged and do not count as support")