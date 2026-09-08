"""
Wrappers to run MAGeCK subcommands as subprocesses.

These functions shell out to the `mageck` binary (installed via
Bioconda in the Docker image) and return structured results instead
of raising raw subprocess exceptions, so the Streamlit UI can display
clear error messages.
"""

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PipelineResult:
    """Outcome of running a MAGeCK subprocess command.

    Attributes:
        success: True if the command exited with code 0.
        stdout: captured standard output.
        stderr: captured standard error (often contains progress logs too).
        output_prefix: the -n prefix used, useful to locate result files.
    """

    success: bool
    stdout: str
    stderr: str
    output_prefix: str


def run_mageck_count(
    fastq_paths: list[str | Path],
    sample_labels: list[str],
    library_path: str | Path,
    output_dir: str | Path,
    output_prefix: str = "mageck_run",
) -> PipelineResult:
    """Run `mageck count` on a set of FASTQ files.

    Args:
        fastq_paths: paths to the FASTQ files, in the same order as sample_labels.
        sample_labels: short names for each sample (e.g. ["control", "treatment"]).
        library_path: path to the sgRNA library CSV file.
        output_dir: directory where output files will be written.
        output_prefix: prefix used for all output file names.

    Returns:
        A PipelineResult with success status and captured logs.
    """
    if len(fastq_paths) != len(sample_labels):
        raise ValueError("fastq_paths and sample_labels must have the same length")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "mageck",
        "count",
        "-l",
        str(library_path),
        "-n",
        str(output_dir / output_prefix),
        "--sample-label",
        ",".join(sample_labels),
        "--fastq",
        *[str(p) for p in fastq_paths],
    ]

    result = subprocess.run(command, capture_output=True, text=True)

    return PipelineResult(
        success=result.returncode == 0,
        stdout=result.stdout,
        stderr=result.stderr,
        output_prefix=output_prefix,
    )

def run_mageck_test(
    count_table_path: str | Path,
    treatment_labels: list[str],
    control_labels: list[str],
    output_dir: str | Path,
    output_prefix: str = "mageck_run",
) -> PipelineResult:
    """Run `mageck test` on a count table produced by `mageck count`.

    Args:
        count_table_path: path to the count.txt file (from run_mageck_count).
        treatment_labels: sample labels to treat as the treatment group.
        control_labels: sample labels to treat as the control group.
        output_dir: directory where output files will be written.
        output_prefix: prefix used for all output file names.

    Returns:
        A PipelineResult with success status and captured logs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "mageck",
        "test",
        "-k",
        str(count_table_path),
        "-t",
        ",".join(treatment_labels),
        "-c",
        ",".join(control_labels),
        "-n",
        str(output_dir / output_prefix),
    ]

    result = subprocess.run(command, capture_output=True, text=True)

    return PipelineResult(
        success=result.returncode == 0,
        stdout=result.stdout,
        stderr=result.stderr,
        output_prefix=output_prefix,
    )

def run_mageck_mle(
    count_table_path: str | Path,
    design_matrix_path: str | Path,
    output_dir: str | Path,
    output_prefix: str = "mageck_run",
) -> PipelineResult:
    """Run `mageck mle` on a count table with a design matrix.

    Args:
        count_table_path: path to the count.txt file (from run_mageck_count).
        design_matrix_path: path to a MAGeCK-formatted design matrix file.
        output_dir: directory where output files will be written.
        output_prefix: prefix used for all output file names.

    Returns:
        A PipelineResult with success status and captured logs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        "mageck",
        "mle",
        "-k",
        str(count_table_path),
        "-d",
        str(design_matrix_path),
        "-n",
        str(output_dir / output_prefix),
    ]

    result = subprocess.run(command, capture_output=True, text=True)

    return PipelineResult(
        success=result.returncode == 0,
        stdout=result.stdout,
        stderr=result.stderr,
        output_prefix=output_prefix,
    )