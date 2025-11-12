# CWI-SSN

## Coverage-Weighted Identity-based Sequence Similarity Network

**CWI-SSN** is a lightweight pipeline for constructing sequence similarity networks (SSNs) from all-vs-all BLAST results.

Instead of relying on EFI-EST’s alignment score, CWI-SSN uses a biologically intuitive metric — the **Coverage-Weighted Identity (CWI)** — defined as: $\text{CWI} = \frac{\text{pident}}{100} \times \sqrt{\frac{L_{align}}{L_Q} \times \frac{L_{align}}{L_S}}$

This measure integrates sequence identity with symmetric coverage, providing a length-independent and easily interpretable similarity score (0–1).

CWI-SSN is particularly suited for **short or highly variable sequences**, such as **RiPP precursor/core peptides**, where traditional statistical alignment scores lose their significance. Therefore, our analysis focuses on sequence identity rather than bit score or E-value.



## 🚀 How to use

**Create environment and install requirements** 

```
conda create -n modifiedSSN python=3.14
conda install biopython blast pandas networkx
```

**Run and build SSN**

```bash
python3 cwi_ssn.py -i task/neuripp_deepripp_trripp_precursor.faa -o task/result --cwi-threshold 0.6 --threads 6 --include-singletons
```



## ❓ Why not EFI-EST SSN for short peptides

The original EFI-EST framework builds sequence similarity networks (SSNs) using **BLAST-based pairwise comparisons**, where edges are scored by an *alignment score* ($AS$) derived from the *bit score* ($B$) ($AS = -\log_{10}\!\left[ 2^{-B} (L_Q L_S) \right]$) . Where, $L_Q$ and $L_S$ are the lengths of the query and subject sequences, respectively.

EFI-EST was originally developed for enzymes of comparable lengths, where a single alignment score threshold reliably correlates with functional similarity. However, for short or length-variable sequences—such as **RiPP precursor or core peptides**—this assumption fails, resulting in biased similarity scores and fragmented networks. The core issue lies in the **length bias of bitscore-derived alignment scores**. According to the definition of *alignment score*, its length dependence originates from the *bit score*, which in turn derives from the *raw alignment score* in BLAST. Because the *raw alignment score* increases proportionally with the length of the aligned region, shorter identical sequences naturally produce lower *bit scores* and *alignment score*, relative to longer sequences. Consequently, even sequences that are compositionally identical but differ in length can yield substantially different *bit* and *alignment scores*. Therefore, a single *alignment scores* threshold is not suitable for short or length-variable sequences.

A small illustrative case is documented in the [issue](issue) file. I have also initiated a discussion regarding this issue at [EFI-EST Discussion #222](https://github.com/EnzymeFunctionInitiative/EST/discussions/222).




