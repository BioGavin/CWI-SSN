# CWI-SSN

## 🎯 Coverage-Weighted Identity-based Sequence Similarity Network

**CWI-SSN** is a lightweight pipeline for constructing sequence similarity networks (SSNs) from all-vs-all BLAST results.

Instead of relying on EFI-EST’s alignment score, CWI-SSN uses a biologically intuitive metric — the **Coverage-Weighted Identity (CWI)** — defined as: $\text{CWI} = \frac{\text{pident}}{100} \times \sqrt{\frac{L_{align}}{L_Q} \times \frac{L_{align}}{L_S}}$

This measure integrates sequence identity with symmetric coverage, providing a length-independent and easily interpretable similarity score (0–1).

CWI-SSN is particularly suited for **short or highly variable sequences**, such as **RiPP precursor/core peptides**, where traditional statistical alignment scores lose their significance. Therefore, our analysis focuses on sequence identity rather than bit score or E-value.



## 🚀 Installation & Quick Start

### Create environment and install requirements
```bash
mamba create -n modifiedSSN python=3.12 blast=2.16.0 -y
mamba activate modifiedSSN
```

### Install from source

```bash
git clone https://github.com/BioGavin/CWI-SSN.git
cd CWI-SSN
pip install .
```

This will install CWI-SSN with all dependencies and create the `cwi_ssn` command-line tool.

**Note:** Always install `biopython` via conda (`mamba install -c bioconda biopython`) rather than pip to avoid compilation issues on some systems.

### Basic Usage

```bash
cwi-ssn -i sequences.faa -o results --cwi-threshold 0.6 --threads 6
```

### Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--input` | `-i` | **required** | Protein FASTA file for BLAST database and queries |
| `--output-dir` | `-o` | `results` | Output directory for BLAST and SSN files |
| `--threads` | | `16` | CPU threads for BLASTP |
| `--cwi-threshold` | | `0.6` | Minimum CWI score (0-1) for keeping edges |
| `--include-singletons` | | `False` | Include sequences without edges as singleton nodes |

### Example

```bash
cwi-ssn -i data/demo.faa -o my_results --cwi-threshold 0.6 --threads 6 --include-singletons
```

### Output Files

The pipeline generates:
- `blast/` - BLAST database and all-vs-all results (`.m8` format)
- `ssn/` - Processed network files:
  - `network_cwi*.graphml` - Final network in GraphML format
  - Intermediate files (reduced edges, scored edges)



## ❓ Why not EFI-EST SSN for short peptides

The original EFI-EST framework builds sequence similarity networks (SSNs) using **BLAST-based pairwise comparisons**, where edges are scored by an *alignment score* ($AS$) derived from the *bit score* ($B$) ($AS = -\log_{10}\!\left[ 2^{-B} (L_Q L_S) \right]$) . Where, $L_Q$ and $L_S$ are the lengths of the query and subject sequences, respectively.

EFI-EST was originally developed for enzymes of comparable lengths, where a single alignment score threshold reliably correlates with functional similarity. However, for short or length-variable sequences—such as **RiPP precursor or core peptides**—this assumption fails, resulting in biased similarity scores and fragmented networks. The core issue lies in the **length bias of bitscore-derived alignment scores**. According to the definition of *alignment score*, its length dependence originates from the *bit score*, which in turn derives from the *raw alignment score* in BLAST. Because the *raw alignment score* increases proportionally with the length of the aligned region, shorter identical sequences naturally produce lower *bit scores* and *alignment score*, relative to longer sequences. Consequently, even sequences that are compositionally identical but differ in length can yield substantially different *bit* and *alignment scores*. Therefore, a single *alignment scores* threshold is not suitable for short or length-variable sequences.

A small illustrative case is documented in the [issue](issue) file. I have also initiated a discussion regarding this issue at [EFI-EST Discussion #222](https://github.com/EnzymeFunctionInitiative/EST/discussions/222).



