import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.fastq_qc import compute_gc_distribution, compute_quality_profile

fixture = Path(__file__).parent / "fixtures_sample.fastq"

gc_df = compute_gc_distribution(fixture)
print("GC distribution:")
print(gc_df)

assert len(gc_df) == 4, f"Expected 4 reads, got {len(gc_df)}"
# READ_1 (GCGCGCGCATATATATGCGC) should have high GC
assert gc_df["gc_percent"].iloc[0] > 50
# READ_4 (AAAATTTTAAAATTTTAAAA) should have 0% GC
assert gc_df["gc_percent"].iloc[3] == 0
print("OK: GC percentages match expected values")

quality_df = compute_quality_profile(fixture)
print("\nQuality profile (first 5 positions):")
print(quality_df.head())

assert len(quality_df) == 20, f"Expected 20 positions, got {len(quality_df)}"
# Position 1 (index 0) should be high quality across most reads ('I' = Q40)
# except READ_4 which is all '!' (Q0) -> average should be pulled down but not to 0
assert 0 < quality_df["mean_quality"].iloc[0] < 40
print("OK: quality profile computed for all positions")

from app.plots import build_gc_content_plot, build_quality_profile_plot

gc_fig = build_gc_content_plot(gc_df, sample_name="fixtures_sample")
assert gc_fig.data, "GC content figure has no traces"
print(f"\nOK: GC content figure built with {len(gc_fig.data)} trace(s)")

quality_fig = build_quality_profile_plot(quality_df, sample_name="fixtures_sample")
assert quality_fig.data, "Quality profile figure has no traces"
print(f"OK: quality profile figure built with {len(quality_fig.data)} trace(s)")