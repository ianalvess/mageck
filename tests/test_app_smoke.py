"""Smoke test for the Streamlit app: every section must render without raising.

Part 1 renders every section with no data loaded. Part 2 injects the small test
fixtures through st.session_state so the "use last run" code paths (plots, metrics,
download buttons) are executed as well. Run it from the project root.
"""

import logging
import sys
from pathlib import Path

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parent.parent
APP_PATH = str(ROOT / "app" / "main.py")
FIXTURES = ROOT / "tests"

# `streamlit run app/main.py` puts the app folder on sys.path, which is why main.py
# uses flat imports (from parser import ...). AppTest does not, so mirror it here.
sys.path.insert(0, str(ROOT / "app"))

# Outside `streamlit run` this logger warns on every run; it is harmless noise.
logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context").setLevel(logging.ERROR)


def assert_no_exception(at: AppTest, label: str) -> None:
    assert not at.exception, f"{label} raised: {[getattr(e, 'value', e) for e in at.exception]}"
    print(f"OK: {label}")


def open_section(section: str, session_state: dict | None = None) -> AppTest:
    at = AppTest.from_file(APP_PATH, default_timeout=120)
    for key, value in (session_state or {}).items():
        at.session_state[key] = value
    at.run()
    assert not at.exception, (
        f"Initial render raised: {[getattr(e, 'value', e) for e in at.exception]}"
    )
    at.sidebar.radio[0].set_value(section).run()
    return at


def use_last_run(at: AppTest, radio_key: str) -> None:
    """Switch a data-source radio to its second option (the "use last ..." choice)."""
    radio = at.radio(key=radio_key)
    radio.set_value(radio.options[1]).run()


# Part 1: every section renders with no data loaded
for section_name in ["Pipeline", "Quality Control", "Results", "Comparisons"]:
    assert_no_exception(open_section(section_name), f"section '{section_name}' renders with no data")

# Part 2: "use last run" paths with the small fixtures
count_table = FIXTURES / "fixtures" / "replicate_counts.txt"

at = open_section(
    "Quality Control",
    {"pipeline_countsummary": FIXTURES / "fixtures" / "countsummary.txt", "pipeline_count_table": count_table},
)
use_last_run(at, "mapping_source")
assert_no_exception(at, "Mapping QC with a countsummary")
use_last_run(at, "qc_source")
assert_no_exception(at, "QC Report with a count table")

at = open_section(
    "Results",
    {
        "pipeline_gene_summary": FIXTURES / "fixtures" / "gene_summary.txt",
        "pipeline_mle_gene_summary": FIXTURES / "fixtures" / "mle_gene_summary.txt",
    },
)
use_last_run(at, "volcano_source")
assert_no_exception(at, "Test Results with a gene_summary")
use_last_run(at, "mle_view_source")
assert_no_exception(at, "MLE Results with an MLE gene_summary")
use_last_run(at, "consensus_rra_source")
use_last_run(at, "consensus_mle_source")
assert_no_exception(at, "Consensus with RRA and MLE results")


at = open_section("Comparisons", {"pipeline_count_table": count_table})
use_last_run(at, "replicates_source")
assert_no_exception(at, "Replicate Correlation with a count table")
use_last_run(at, "convergence_source")
at.multiselect(key="convergence_controls").set_value(["control_rep1", "control_rep2"]).run()
at.multiselect(key="convergence_replicates").set_value(["treatment_rep1", "treatment_rep2"]).run()
assert_no_exception(at, "Replicate Convergence with a count table")

print("OK: all sections render without exceptions")