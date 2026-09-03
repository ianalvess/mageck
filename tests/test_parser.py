import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_gene_summary

fixture = Path(__file__).parent / "fixtures_gene_summary.txt"
df = load_gene_summary(fixture)

print(df[["id", "direction", "best_fdr", "best_lfc"]])

assert df.loc[df["id"] == "GENE_A", "direction"].iloc[0] == "neg"
assert df.loc[df["id"] == "GENE_B", "direction"].iloc[0] == "pos"
print("\nOK: direction escolhida corretamente para GENE_A (neg) e GENE_B (pos)")