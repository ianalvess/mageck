import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from app.plots import build_correlation_heatmap, build_pca_plot

fixture = Path(__file__).parent / "fixtures_replicate_counts.txt"
df = pd.read_csv(fixture, sep="\t")

sample_columns = ["control_rep1", "control_rep2", "treatment_rep1", "treatment_rep2"]

heatmap_fig = build_correlation_heatmap(df, sample_columns)
assert heatmap_fig.data, "Heatmap figure has no traces"
print(f"OK: correlation heatmap built with {len(heatmap_fig.data)} trace(s)")

pca_fig = build_pca_plot(df, sample_columns)
assert pca_fig.data, "PCA figure has no traces"
print(f"OK: PCA figure built with {len(pca_fig.data)} trace(s)")

# Sanity check: control replicates should be closer to each other in PC space
# than to treatment replicates, given how the fixture was designed
pca_data = pca_fig.data[0]
print("PCA coordinates:")
for name, x, y in zip(pca_data.text, pca_data.x, pca_data.y):
    print(f"  {name}: PC1={x:.2f}, PC2={y:.2f}")