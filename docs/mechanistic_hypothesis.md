# Mechanistic Hypothesis: Immune Regulatory Balance in Symptomatic vs Asymptomatic Malaria

## 1. Summary of Problem and Data

Malaria caused by *Plasmodium falciparum* can produce vastly different clinical outcomes: some infected individuals develop fever and severe symptoms, while others carry the parasite without noticeable illness. Understanding what drives these divergent states is critical for vaccine design and elimination efforts.

We analyzed public RNA‑seq data from PBMCs (GSE181179) comparing symptomatic (n=6), asymptomatic (n=5), and healthy controls (n=6). Using limma‑voom with global FDR correction (as in the original study by Studniberg et al., 2022), we replicated the original differential expression results almost exactly: 928 DEGs for Symptomatic vs Healthy, 178 for Asymptomatic vs Healthy, and 417 for the key comparison Symptomatic vs Asymptomatic.

This document presents a mechanistic hypothesis based on those results, aiming to explain how the same parasite can lead to either clinical malaria or a silent, regulated infection.

## 2. Key Findings from Our Analysis

Our pathway and network analyses revealed a clear dichotomy:

- **Up‑regulated in symptomatic individuals** (relative to asymptomatic): genes involved in cell cycle, mitosis, DNA replication, and inflammatory signaling. Top hubs include PSMC2, PSMD14, E2F1, PCNA, and multiple histone genes (H2BC, H2AC). These reflect a proliferative and stress‑like state in immune cells.

- **Down‑regulated in symptomatic individuals (i.e., up in asymptomatic)**: genes with known anti‑proliferative or immunoregulatory functions, including SIRT1, SMAD7, BTG1, IRF2BP2, and NR1D1. These are not simply passive markers; they actively suppress inflammation and T‑cell activation.

The top gene hubs by network degree:

| Gene   | Degree | logFC (SM vs AM) |
|--------|--------|------------------|
| SIRT1  | 253    | -0.61            |
| CAV1   | 208    | 1.69             |
| PSMD14 | 169    | 0.52             |
| PSMC2  | 162    | 0.55             |
| EGF    | 158    | 3.63             |
| PSMA5  | 156    | 0.53             |
| CLU    | 133    | 2.97             |
| E2F1   | 116    | 1.98             |
| H2BC11 | 115    | 2.48             |
| ITGB3  | 108    | 2.58             |

While the up‑regulated hubs are dominated by proliferative genes, the presence of SIRT1 and SMAD7 as down‑regulated hubs (i.e., higher in asymptomatic) is notable. SIRT1 deacetylates NF‑κB and p53, reducing inflammatory gene expression. SMAD7 inhibits TGF‑β and NF‑κB signaling, acting as a brake on immune activation.

## 3. Proposed Molecular Mechanism

We propose that the clinical state is determined by the net balance between two opposing transcriptional modules:

1. **Pro‑inflammatory / proliferative module**  
   - Parasite products (e.g., hemozoin, GPI anchors) are sensed by innate receptors (TLRs).  
   - This activates NF‑κB and MAPK signaling cascades.  
   - Downstream cytokines (TNF, IL‑6, IFN‑γ) and cell cycle regulators (E2F1, PCNA, cyclins) drive immune cell activation and proliferation.  

2. **Regulatory / anti‑proliferative module**  
   - Immune checkpoints (CTLA‑4) and negative regulators (SIRT1, SMAD7, BTG1, IRF2BP2) suppress T‑cell activation and cell cycle progression.  
   - These factors are upregulated in asymptomatic malaria, damping inflammation and allowing parasite persistence without clinical symptoms.

A simple conceptual diagram:

    Malaria PAMPs
        ↓
    Innate sensing (TLRs)
        ↓
    ┌───────────────┐
    │ NF-κB / MAPK │ → Inflammatory cytokines (TNF, IL-6, IFN-γ)
    └───────────────┘
        ↓
    T-cell activation & proliferation (cell cycle genes)
        ↓
    Regulatory feedback (CTLA-4, SIRT1, SMAD7, IL-10)
        ↓
    Suppression of proliferation / inflammation

In symptomatic infection, the pro‑inflammatory module dominates, leading to fever, immunopathology, and high parasitemia. In asymptomatic infection, the regulatory module is upregulated, suppressing the inflammatory response and limiting clinical symptoms, at the cost of incomplete parasite clearance.

## 4. Key Molecular Players for the Model

For a minimal dynamic model, we will represent the system with the following nodes:

- **NF‑κB** – pro‑inflammatory transcription factor  
- **MAPK** – stress/cytokine signaling kinase  
- **Inflammatory output** – represented by cytokines (e.g., IFN‑γ, TNF)  
- **CTLA‑4** – immune checkpoint receptor, inhibits T‑cell activation  
- **SIRT1** – deacetylase, negative regulator of NF‑κB  
- **SMAD7** – inhibitor of TGF‑β/NF‑κB signaling  
- **IL‑10** (optional) – anti‑inflammatory cytokine, can be included later

The model will explore how changes in these regulators shift the system from a regulated (asymptomatic) to an inflammatory (symptomatic) state.

## 5. Assumptions and Simplifications

We acknowledge several important limitations:

- Bulk RNA‑seq provides an average snapshot of gene expression across mixed PBMC populations. It cannot resolve cell‑type‑specific dynamics or post‑translational modifications.
- The causal relationships proposed here are based on existing literature, not directly proven by our transcriptomic data alone.
- We assume that differences in mRNA levels correspond to functional differences at the protein level.
- The model will be **qualitative or semi‑quantitative**: parameters will be chosen to reproduce the observed qualitative behavior, not fitted to quantitative data.

Therefore, this hypothesis is a **starting point for model development and experimental testing**, not a definitive explanation.

## 6. Testable Predictions

If the proposed mechanism is correct, the model should make the following predictions:

1. **Enhancing CTLA‑4 or SIRT1 activity** should shift the system toward an asymptomatic‑like state (reduced proliferation/inflammation).  
2. **Inhibiting NF‑κB or MAPK** should reduce inflammatory output and cell cycle gene expression.  
3. **Perturbing SMAD7** should alter the balance between TGF‑β and NF‑κB signaling, affecting T‑cell proliferation.  
4. The model should reproduce the qualitative gene expression patterns: high inflammatory/proliferative signature in symptomatic, high regulatory signature in asymptomatic.

These predictions can be tested experimentally in PBMC cultures or mouse models (e.g., CTLA‑4 blockade, SIRT1 agonists).

## 7. Next Steps: From Hypothesis to Model

The next phase will translate this conceptual model into a system of ordinary differential equations (ODEs). We will:

- Define state variables for NF‑κB, MAPK, inflammatory output, and regulatory feedback.
- Assign biologically plausible parameters (based on literature estimates or reasonable ranges).
- Simulate baseline and perturbed conditions (e.g., increased CTLA‑4, reduced NF‑κB).
- Perform sensitivity analysis to identify critical nodes.

This will allow us to test whether a simple regulatory network can account for the two clinical states observed in malaria. The model will remain openly documented and reproducible, as part of this project.