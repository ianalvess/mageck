import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_mle_gene_summary
from app.plots import build_beta_comparison_plot

fixture = Path(__file__).parent / "fixtures_mle_gene_summary.txt"
df, conditions = load_mle_gene_summary(fixture)

print("Conditions found:", conditions)
print(df[["Gene", "PLX|beta", "PLX|fdr"]])

assert conditions == ["Control", "PLX"]
assert len(df) == 4
print("\nOK: conditions and row count match expected values")

fig = build_beta_comparison_plot(df, condition_x="Control", condition_y="PLX", fdr_threshold=0.05)
assert fig.data, "Figure has no traces"
print(f"\nOK: beta comparison figure built with {len(fig.data)} trace(s)")

sig_both = [pt for trace in fig.data if trace.name == "Both" for pt in trace.hovertext]
print("Significant in both:", sig_both)
assert "GENE_D" in sig_both
print("OK: GENE_D correctly flagged as significant in both conditions")