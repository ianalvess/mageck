"""
Lightweight FASTQ quality-control utilities.

These functions stream reads from a (potentially very large) FASTQ
file to compute GC content and per-base quality statistics, without
loading the whole file into memory at once.
"""

import array
import gzip
import statistics
from pathlib import Path

import pandas as pd


def _open_fastq(path: str | Path):
    """Open a FASTQ file, transparently handling gzip compression."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def sample_fastq_reads(path: str | Path, max_reads: int | None = None):
    """Yield (sequence, quality_string) tuples for reads in a FASTQ file.

    Args:
        path: path to a FASTQ (optionally gzip-compressed) file.
        max_reads: if set, stop after this many reads. If None (default),
            process the entire file.

    Yields:
        Tuples of (sequence, quality_string) for each read.
    """
    with _open_fastq(path) as f:
        count = 0
        while max_reads is None or count < max_reads:
            header = f.readline()
            if not header:
                break  # end of file reached
            sequence = f.readline().strip()
            f.readline()
            quality = f.readline().strip()
            yield sequence, quality
            count += 1


def compute_gc_distribution(path: str | Path, max_reads: int | None = None) -> pd.DataFrame:
    """Compute per-read GC content across a FASTQ file.

    Args:
        path: path to a FASTQ file.
        max_reads: if set, only process this many reads; otherwise process all.

    Returns:
        A DataFrame with one row per read and a "gc_percent" column.
    """
    gc_percentages = []
    for sequence, _quality in sample_fastq_reads(path, max_reads=max_reads):
        if not sequence:
            continue
        gc_count = sequence.upper().count("G") + sequence.upper().count("C")
        gc_percentages.append(100 * gc_count / len(sequence))

    return pd.DataFrame({"gc_percent": gc_percentages})


def compute_quality_profile(path: str | Path, max_reads: int | None = None) -> pd.DataFrame:
    """Compute median Phred quality score at each base position, across a FASTQ file.

    Assumes Phred+33 encoding (the standard for modern Illumina data).
    Uses a byte array per position (1 byte per quality score) to keep
    memory usage as low as possible for very large files, since the
    true median requires keeping every value in memory (unlike a
    running mean).

    Args:
        path: path to a FASTQ file.
        max_reads: if set, only process this many reads; otherwise process all.

    Returns:
        A DataFrame with columns "position" (1-indexed) and "median_quality".
    """
    position_values: dict[int, array.array] = {}

    for _sequence, quality in sample_fastq_reads(path, max_reads=max_reads):
        for i, char in enumerate(quality):
            phred_score = ord(char) - 33
            if i not in position_values:
                position_values[i] = array.array("B")
            position_values[i].append(phred_score)

    positions = sorted(position_values.keys())
    median_qualities = [statistics.median(position_values[p]) for p in positions]

    return pd.DataFrame({"position": [p + 1 for p in positions], "median_quality": median_qualities})