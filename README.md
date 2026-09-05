# Malaria Host-Response Systems Explorer

**A computational investigation into the transcriptional and regulatory differences between symptomatic and asymptomatic *Plasmodium falciparum* malaria.**

## Overview

This project uses public transcriptomic data (GSE181179) to uncover the molecular mechanisms that distinguish symptomatic from asymptomatic malaria. By integrating differential expression, pathway enrichment, network analysis, and eventually dynamic modeling, we aim to identify key regulatory hubs that control the balance between inflammation and immune suppression.

The overarching question is:

> How do transcriptional states of human immune cells differ between symptomatic malaria, asymptomatic malaria, and healthy controls, and what biological pathways and regulatory relationships may help explain those differences?

## Current Status

- **Phase 1–3**: Biological question map complete (`docs/biological_question.md`).
- **Phase 4–6**: Raw data merged, metadata created, and count matrix validated.
- **Phase 7**: Quality control and exploratory analysis (PCA, sample distance, heatmap).
- **Phase 8–9**: Differential expression with limma-voom; results replicate original study (Studniberg et al., *Mol Syst Biol*, 2022).
- **Phase 10**: Pathway enrichment (GO, KEGG, Reactome) completed.
- **Phase 11–12**: Gene–pathway network analysis; identified candidate regulatory hubs.
- **Phase 13**: Mechanistic hypothesis formulated (`docs/mechanistic_hypothesis.md`).
- **Next**: Dynamic ODE model of the NF-κB/MAPK/immune checkpoint axis.

## Repository Structure
malaria-host-response-systems/
├── data/
│ ├── raw/ # Downloaded per-sample count files (tar, extracted)
│ └── processed/ # Merged counts, metadata, differential expression results
├── notebooks/
│ ├── 01_data_exploration.ipynb
│ ├── 02_differential_expression.ipynb # (legacy, PyDESeq2 attempt)
│ ├── 03_gene_annotation_and_plots.ipynb # PCA, sample distance, heatmap
│ ├── 04_pathway_analysis.ipynb # Enrichment analysis
│ └── 05_network_analysis.ipynb # Gene–pathway network
├── src/ # Reusable Python modules (future)
├── models/ # Dynamic models (future)
├── results/
│ └── figures/ # Generated plots (PNG/PDF)
├── docs/ # Biological question, mechanistic hypothesis
├── run_limma_voom.R # R script for differential expression (limma-voom)
├── merge_counts.py # Python script to merge per-sample count files
└── README.md


## Key Findings So Far

- **Symptomatic malaria** is associated with up-regulation of cell cycle, DNA replication, and inflammatory pathways.
- **Asymptomatic malaria** shows up-regulation of negative regulators of immune activation, including *SIRT1*, *SMAD7*, *BTG1*, and *IRF2BP2*.
- The balance between pro-inflammatory/proliferative signals (NF-κB, MAPK) and regulatory signals (CTLA-4, SIRT1, SMAD7) may determine clinical outcome.

## Reproducing the Analysis

### 1. Environment

The analysis uses both R and Python. Key R packages: `limma`, `edgeR`, `DESeq2` (optional). Python packages: `pandas`, `numpy`, `matplotlib`, `seaborn`, `gseapy`, `networkx`, `pyvis`.

### 2. Data Acquisition

Download the raw count files from GEO (GSE181179). Extract the tar archive and merge individual sample files using `merge_counts.py`.

### 3. Differential Expression
Run the R script:
```bash Rscript run_limma_voom.R ```bash

This produces CSV files of DE results for all three comparisons.

### 4. Downstream Analysis
Use the notebooks in order (03_*, 04_*, 05_*) to perform QC, pathway enrichment, and network analysis.

Authors
Alule Robert
MSc Biochemistry student, Makerere University
Research Interests: computational biology, immunology, host–pathogen interactions, systems modeling

License
This project is licensed under the MIT License.





