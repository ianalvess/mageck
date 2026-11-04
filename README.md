# MAGeCK Dashboard



A Dockerized Streamlit dashboard that makes CRISPR screen analysis with \[MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) easier. Run `mageck count`, `test` and `mle` from the browser, check data quality, and explore the results with interactive plots.



## Features



**Pipeline**: run `mageck count`, `mageck test` and `mageck mle` on files that live in `data/input/`, and compute GC content and base quality directly from raw FASTQ files.

**Quality Control**: mapping rate and zero-count sgRNAs per sample, per-sample read distribution with two Gini indices, and a sample overview (correlation heatmap and PCA).

**Results**: volcano plot for `mageck test`, beta scores for `mageck mle` (permutation or Wald FDR), and a consensus table that counts how many methods (RRA, MLE, replicate convergence) support each gene.

**Comparisons**: replicate correlation, comparison of two screens, and replicate convergence (sgRNAs enriched in every replicate).

**Helpers**: detection of sgRNAs that were never observed (e.g. an unused half of a split library), validation of count tables, and PNG export of every chart.



## Quick start



```bash

git clone https://github.com/ianalvess/mageck.git

cd mageck

docker compose up --build

```



Open http://localhost:8501. The `data/` folder is mounted into the container.



## Typical workflow



1\. Copy your FASTQ files and the sgRNA library CSV into `data/input/` (subfolders are fine).

2\. **Pipeline > Run Count**: choose the library and one FASTQ per sample, name the samples, and run `mageck count`.

3\. **Quality Control**: check Mapping QC, QC Report and Sample Overview before interpreting any hit.

4\. **Pipeline > Run Test** (two groups) or **Run MLE** (design matrix), then open **Results**.

5\. Use **Comparisons** to check replicate reproducibility and to look for convergent sgRNAs.



Outputs of runs started from the dashboard are written to `data/pipeline\_runs/<prefix>/`.



## Input formats



**Library**: CSV without header: sgRNA ID, sequence, gene.

**Count table**: the `\*.count.txt` from `mageck count` (`sgRNA`, `Gene`, then one numeric column per sample).

**Design matrix** (MLE): tab-separated, header `Samples  baseline  <condition> ...`; the first column must match the sample labels of the count table.



## Interpretation notes



**Split libraries**: sgRNAs with zero reads in every sample usually mean part of the library file was not used. They inflate zero-count statistics, so the dashboard warns about them and can exclude them from the QC and from `test`/`mle` runs.

**Two Gini indices**: the "log2" value follows the MAGeCK countsummary definition (computed on `log2(count + 1)`); the "raw" value is dominated by a few very abundant sgRNAs and flags clonal expansion.

**MLE FDR**: permutation FDR can have a floor (no gene can pass the threshold) when no control sgRNAs are given; Wald FDR can be over-optimistic with few replicates. Compare both, and confirm hits with another method.

**Control sgRNAs**: normalizing with non-targeting controls can change the number of hits dramatically, so treat the extra hits as exploratory.

**Clonal tumors**: when one replicate is dominated by a few clones, gene-level tests can be driven by a single sample. The convergence analysis looks for sgRNAs enriched in every replicate instead.

**Consensus**: ranks candidates by how many methods support them; it does not validate them, because all methods use the same replicates. MLE betas above 10 in absolute value are flagged as unstable and do not count as support.



## Development



The `mageck` binary only exists inside the Docker image, so runs started from the dashboard need Docker. Everything else can run locally:



```bash

pip install -r requirements.txt

streamlit run app/main.py

```



Each file in `tests/` is a plain script. Run them from the project root, e.g. `python tests/test\_parser.py`. `python tests/test\_app\_smoke.py` renders every section of the app headlessly. `python -m pyflakes app tests` catches undefined names.



`numpy` is pinned to 1.26.4 for MAGeCK compatibility, so `pyarrow` is pinned to a release that still supports NumPy 1.x.



## Project layout



```

project/
│
├── app/
│   ├── main.py             # Streamlit interface (four main sections)
│   ├── parser.py           # Readers and parsers for MAGeCK outputs
│   ├── plots.py            # Plotly visualizations
│   ├── pipeline.py         # Wrappers for MAGeCK count, test and MLE
│   ├── fastq_qc.py         # Raw FASTQ quality control (GC content, base quality)
│   ├── library_check.py    # sgRNA library validation and never-observed sgRNA detection
│   ├── convergence.py      # Replicate convergence analysis
│   └── consensus.py        # Gene-level consensus across RRA, MLE and convergence
│
├── tests/
│   └── fixtures/           # Small input files for testing
│
├── data/
│   ├── input/              # FASTQ files, sgRNA libraries and design matrices
│   └── pipeline_runs/      # Outputs generated by dashboard pipeline runs
│
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── .streamlit/
    └── config.toml

```



## License



See [LICENSE](LICENSE).

