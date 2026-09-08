import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_mle_gene_summary
from app.plots import build_beta_score_plot

fixture = Path(__file__).parent / "fixtures_mle_gene_summary.txt"
df, conditions = load_mle_gene_summary(fixture)

print("Conditions found:", conditions)
print(df[["Gene", "PLX|beta", "PLX|fdr"]])

assert conditions == ["PLX"]
assert len(df) == 4
print("\nOK: conditions and row count match expected values")

fig = build_beta_score_plot(df, condition="PLX", fdr_threshold=0.05)
assert fig.data, "Figure has no traces"
print(f"\nOK: beta score figure built with {len(fig.data)} trace(s)")

sig_genes = [pt for trace in fig.data if trace.name == "Significant" for pt in trace.hovertext]
print("Significant genes:", sig_genes)
assert "GENE_A" in sig_genes
assert "GENE_B" in sig_genes
assert "GENE_C" not in sig_genes
print("OK: significance flags match expected genes")