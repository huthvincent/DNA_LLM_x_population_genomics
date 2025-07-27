# 🧬 Harnessing DNA Foundation Models for Human Population Genomics  
### Insights from the **APOE** Locus   |   Science submission code

> **One-liner**  
> **Evo2** DNA foundation model + population-scale WGS → variant & haplotype scores that illuminate Alzheimer’s genetics across ADNI, HPRC & UK Biobank.

---

## 📖 Table of Contents
1. [Project Snapshot](#project-snapshot)
2. [Quick Start](#quick-start)
3. [Repository Layout](#repository-layout)
4. [Reproduce the Figures](#reproduce-the-figures)
5. [Data & Model](#data--model)
6. [License](#license)
7. [Contact](#contact)

---

## Project Snapshot
DNA foundation models (e.g. **Evo2-7B**) learn regulatory grammar directly from sequence.  
Here we score every variant *and* full-length haplotype in the APOE region (chr19:44.84–44.92 Mb, GRCh38) and show:

| Figure | Highlight | Dataset |
| :----: | --------- | ------- |
| **2** | Δ-scores track AD GWAS peaks | IGAP / ADNI GWAS |
| **3** | Haplotype diversity & ancestry effects | Human Pangenome |
| **4** | Scores predict amyloid & cognition | ADNI WGS + PET |
| **5** | MRI, proteomic & PheWAS landscape | UK Biobank |

The full manuscript title is **“Harnessing DNA Foundation Models for Human Population Genomics: Insights from the APOE Locus.”**

---

## Quick Start
```bash
# 1️⃣ Clone + install R deps (renv for reproducibility)
git clone https://github.com/your-org/APOE-Evo2.git
cd APOE-Evo2
Rscript -e "install.packages('renv'); renv::restore()"

# 2️⃣ Point to data & Evo2 model weight directory
export APOE_DATA_DIR=/path/to/inputs
export EVO2_MODEL_DIR=/path/to/evo2-7b

# 3️⃣ Re-run Figure 2 (as example)
Rscript Figure2/0.EVO2_scoring_GWAS.R
Rscript Figure2/1.Fig2A.R
