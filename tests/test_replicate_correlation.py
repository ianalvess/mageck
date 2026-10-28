import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from app.plots import build_replicate_correlation_plot

fixture = Path(__file__).parent / "fixtures_replicate_counts.txt"
df = pd.read_csv(fixture, sep="\t")

fig = build_replicate_correlation_plot(df, "control_rep1", "control_rep2")

assert fig.data, "Figure has no traces"
print(f"OK: figure built with {len(fig.data)} trace(s)")
print("Title:", fig.layout.title.text)

correlation = df["control_rep1"].corr(df["control_rep2"])
assert correlation > 0.9, f"Expected high correlation for replicates, got {correlation:.3f}"
print(f"OK: correlation is high as expected ({correlation:.3f})")

real_format = df.rename(columns={"sgrna": "sgRNA"})
real_fig = build_replicate_correlation_plot(real_format, "control_rep1", "control_rep2")
assert list(real_fig.data[0].hovertext)[:2] == ["sgRNA_1", "sgRNA_2"]
print("OK: hover shows sgRNA IDs for MAGeCK-style 'sgRNA' columns")