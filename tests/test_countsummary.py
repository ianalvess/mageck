import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.parser import load_countsummary

from app.plots import build_zerocounts_plot

from app.plots import build_mapping_summary_plot

fixture = Path(__file__).parent / "fixtures" / "countsummary.txt"
df = load_countsummary(fixture)

print(df[["Label", "Percentage", "Zerocounts", "GiniIndex"]])

assert len(df) == 4, f"Expected 4 samples, got {len(df)}"
assert df["Percentage"].between(0, 1).all(), "Percentage values out of [0, 1] range"

low_mapping = df[df["Percentage"] < 0.8]
print(f"\nSamples with mapping < 80%: {low_mapping['Label'].tolist()}")
assert set(low_mapping["Label"]) == {"treatment_rep1", "treatment_rep2"}
print("OK: low-mapping samples correctly identified")

zerocounts_fig = build_zerocounts_plot(df)
assert zerocounts_fig.data, "Figure has no traces"
print(f"OK: zerocounts figure built with {len(zerocounts_fig.data)} trace(s)")

mapping_fig = build_mapping_summary_plot(df, threshold=0.6)
assert any("60%" in annotation.text for annotation in mapping_fig.layout.annotations)
print("OK: mapping plot threshold line follows the threshold argument")