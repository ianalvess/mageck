import sys
from pathlib import Path
import numpy as np
import pandas as pd

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

gini_log = compute_gini_index(df["control_count"], log_scale=True)
assert 0 <= gini_log <= 1, "Log-scale Gini index out of expected [0, 1] range"
assert np.isnan(compute_gini_index(pd.Series([0, 0, 0]))), "All-zero sample should give NaN"
print(f"OK: log-scale Gini = {gini_log:.3f}, all-zero sample gives NaN")
