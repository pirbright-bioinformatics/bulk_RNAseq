# RNAseq QC and Quantification Pipeline

This pipeline provides automated Quality Control (QC), trimming, and quantification for RNAseq data. It supports both paired-end and single-end data and automatically replicates your input directory structure in the output.

## Features

- **Zero-Config Discovery:** Automatically finds samples by scanning your data directory. No sample sheet required.
- **Path-Based Metadata:** Your directory structure (e.g., `condition/tissue/time/animal/`) is preserved and mirrored in the results.
- **SE/PE Support:** Automatically detects if a sample is Single-End or Paired-End based on the number of FASTQ files in its directory.
- **Tools:** FastQC, MultiQC, fastp, and Kallisto.

## Components

1.  **Snakemake.qc**: Handles Quality Control and Trimming.
    - Runs FastQC on raw and trimmed data.
    - Performs trimming using `fastp`.
    - Generates MultiQC reports.
2.  **Snakefile.kallisto**: Handles quantification using Kallisto.
    - Quantifies transcript abundance.
    - Generates abundance plots.

## Usage

### 1. Configuration

The pipeline is configured via `config.yaml`.

```yaml
sampledir: path/to/your/raw_data  # Root directory of your FASTQ files
outdir: results/                 # Where results will be saved
qctype: fastp                    # Trimming method
adapters: adapters.fa            # Path to adapters for fastp

kallisto:
  index: index/index.idx         # Kallisto index path
  # Required for Single-End data:
  fragment_len: 200              
  sd: 20                         
```

### 2. Input Structure

Organize your data into directories. Each **leaf directory** containing FASTQ files (`.fq.gz` or `.fastq.gz`) is treated as one sample.

Example:
```
data/
└── mock/
    └── throat/
        └── 24h/
            └── animal_001/
                ├── read_1.fq.gz
                └── read_2.fq.gz  # 2 files = Paired-End
```

### 3. Execution

Run QC and Trimming:
```bash
snakemake -s Snakemake.qc --configfile config.yaml -j 8
```

Run Quantification:
```bash
snakemake -s Snakefile.kallisto --configfile config.yaml -j 8
```

## Output Structure

The results mirror your input hierarchy:
- `results/fastqc/mock/throat/24h/animal_001/`
- `results/fastp/mock/throat/24h/animal_001/`
- `results/kallisto/mock/throat/24h/animal_001/`
