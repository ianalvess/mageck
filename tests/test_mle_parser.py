import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_mle_gene_summary
from app.plots import build_beta_comparison_plot
from app.plots import build_beta_score_plot

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

single_fig = build_beta_score_plot(df, condition="PLX", fdr_threshold=0.05)
sig_genes = [pt for trace in single_fig.data if trace.name == "Significant" for pt in trace.hovertext]
assert set(sig_genes) == {"GENE_A", "GENE_B", "GENE_D"}
print("OK: single-condition beta score plot flags the expected genes")

zero_df = df.copy()
zero_df.loc[0, "PLX|fdr"] = 0.0
zero_fig = build_beta_score_plot(zero_df, condition="PLX")
max_y = max(max(trace.y) for trace in zero_fig.data)
assert max_y < 10, f"FDR=0 should be floored, got y={max_y}"
print("OK: FDR of 0 no longer produces an extreme y value")

wald_df = df.copy()
wald_df["PLX|wald-fdr"] = [0.5, 0.5, 0.5, 0.001]
wald_fig = build_beta_score_plot(wald_df, condition="PLX", fdr_suffix="wald-fdr")
wald_sig = [pt for trace in wald_fig.data if trace.name == "Significant" for pt in trace.hovertext]
assert list(wald_sig) == ["GENE_D"]
print("OK: fdr_suffix selects the Wald FDR column")