import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_gene_summary
from app.plots import compare_gene_summaries, build_screen_comparison_plot

fixture1 = Path(__file__).parent / "fixtures" / "gene_summary.txt"
fixture2 = Path(__file__).parent / "fixtures" / "gene_summary_screen2.txt"

df1 = load_gene_summary(fixture1)
df2 = load_gene_summary(fixture2)

merged = compare_gene_summaries(df1, df2, fdr_threshold=0.05, label1="Screen1", label2="Screen2")
print(merged[["id", "concordance"]])

# GENE_A should be significant in both -> "Both"
assert merged.loc[merged["id"] == "GENE_A", "concordance"].iloc[0] == "Both"
# GENE_B was significant only in screen 1 -> "Screen1 only"
assert merged.loc[merged["id"] == "GENE_B", "concordance"].iloc[0] == "Screen1 only"
# GENE_D was significant only in screen 2 -> "Screen2 only"
assert merged.loc[merged["id"] == "GENE_D", "concordance"].iloc[0] == "Screen2 only"
print("\nOK: concordance categories match expected values")

fig = build_screen_comparison_plot(merged, label1="Screen1", label2="Screen2")
assert fig.data, "Figure has no traces"
print(f"OK: figure built with {len(fig.data)} trace(s)")