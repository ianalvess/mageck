import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.pipeline import run_mageck_count, run_mageck_test, run_mageck_mle

# Fake a successful subprocess result, without actually running mageck
fake_result = MagicMock()
fake_result.returncode = 0
fake_result.stdout = "done"
fake_result.stderr = ""

with patch("app.pipeline.subprocess.run", return_value=fake_result) as mock_run:
    result = run_mageck_count(
        fastq_paths=["control.fastq", "treatment.fastq"],
        sample_labels=["control", "treatment"],
        library_path="library.csv",
        output_dir="tests/tmp_output",
        output_prefix="test_run",
    )

    assert result.success is True
    assert result.output_prefix == "test_run"

    called_command = mock_run.call_args[0][0]
    print("Command built:", " ".join(called_command))

    assert called_command[0] == "mageck"
    assert called_command[1] == "count"
    assert "--sample-label" in called_command
    assert "control,treatment" in called_command
    assert "control.fastq" in called_command
    assert "treatment.fastq" in called_command

print("\nOK: command constructed correctly and result parsed correctly")

# Test the mismatched-length validation
try:
    run_mageck_count(
        fastq_paths=["only_one.fastq"],
        sample_labels=["control", "treatment"],
        library_path="library.csv",
        output_dir="tests/tmp_output",
    )
    raise AssertionError("Expected ValueError for mismatched lengths")
except ValueError:
    print("OK: mismatched fastq_paths/sample_labels correctly raises ValueError")

# Test run_mageck_test with a mocked subprocess call
with patch("app.pipeline.subprocess.run", return_value=fake_result) as mock_run:
    result = run_mageck_test(
        count_table_path="tests/tmp_output/test_run.count.txt",
        treatment_labels=["treatment"],
        control_labels=["control"],
        output_dir="tests/tmp_output",
        output_prefix="test_run_analysis",
    )

    assert result.success is True

    called_command = mock_run.call_args[0][0]
    print("\nCommand built:", " ".join(called_command))

    assert called_command[0] == "mageck"
    assert called_command[1] == "test"
    assert "-t" in called_command
    assert "treatment" in called_command
    assert "-c" in called_command
    assert "control" in called_command

print("OK: mageck test command constructed correctly")

# Test run_mageck_mle with a mocked subprocess call
with patch("app.pipeline.subprocess.run", return_value=fake_result) as mock_run:
    result = run_mageck_mle(
        count_table_path="tests/tmp_output/test_run.count.txt",
        design_matrix_path="tests/tmp_output/design.txt",
        output_dir="tests/tmp_output",
        output_prefix="test_run_mle",
    )

    assert result.success is True

    called_command = mock_run.call_args[0][0]
    print("\nCommand built:", " ".join(called_command))

    assert called_command[0] == "mageck"
    assert called_command[1] == "mle"
    assert "-d" in called_command
    assert "design.txt" in called_command[called_command.index("-d") + 1]

print("OK: mageck mle command constructed correctly")