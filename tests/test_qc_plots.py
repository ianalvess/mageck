import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_sgrna_summary
from app.plots import compute_gini_index, build_read_distribution_plot

fixture = Path(__file__).parent / "fixtures_sgrna_summary.txt"
df = load_sgrna_summary(fixture)

gini = compute_gini_index(df["control_count"])
print(f"Gini index (control_count): {gini:.3f}")
assert 0 <= gini <= 1, "Gini index out of expected [0, 1] range"

fig = build_read_distribution_plot(df, "control_count")
assert fig.data, "Figure has no traces"
print(f"OK: histogram built with {len(fig.data)} trace(s)")