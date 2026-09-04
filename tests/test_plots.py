import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_gene_summary
from app.plots import build_volcano_plot

fixture = Path(__file__).parent / "fixtures_gene_summary.txt"
df = load_gene_summary(fixture)

fig = build_volcano_plot(df, fdr_threshold=0.05, lfc_threshold=1.0)

assert fig.data, "Figure has no traces"
print(f"OK: figure built with {len(fig.data)} trace(s)")

# GENE_A and GENE_B should be flagged significant; GENE_C should not
sig_points = [pt for trace in fig.data if trace.name == "Significant" for pt in trace.hovertext]
print("Significant genes:", sig_points)
assert "GENE_A" in sig_points
assert "GENE_B" in sig_points
assert "GENE_C" not in sig_points
print("OK: significance flags match expected genes")