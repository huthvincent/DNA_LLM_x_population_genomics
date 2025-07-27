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
```
Tip: Each folder’s 0* script builds inputs (variant ↦ FASTA ↦ Δ-score).
Subsequent scripts generate statistics & plots.

⸻

Repository Layout

Path	Role
Figure2/	Variant-level scoring & AD GWAS overlay
Figure3/	Pangenome haplotypes • ancestry correlations
Figure4/	ADNI haplotypes → amyloid & cognition
Figure5/	UKB scale associations (MRI, proteomics, PheWAS)
utils/(optional)	helper scripts for Evo2 API, plotting themes


⸻

## 🔁 Reproduce the Figures

> **One folder = one figure.**  
> Each folder’s **`0*` / `00*`** script builds input sequences & Evo2 Δ-scores.  
> Subsequent scripts analyse and plot. 👇

| Figure | Folder | Quick-run Commands* | Main Outputs |
| :---: | :--- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :---------------- |
| **2** | `Figure2/` | ```bash<br># Variant scoring & GWAS overlay<br>cd Figure2<br>Rscript 0.EVO2_scoring_GWAS.R      # SNP → FASTA → Δ-score<br>Rscript 1.Fig2A.R              # regional plot<br>Rscript 2.Fig2B.R 3.Fig2C.R   # score vs GWAS<br>``` | Manhattan & scatter plots |
| **3** | `Figure3/` | ```bash<br># Pangenome haplotypes<br>bash 00.EVO2_scoring_PanGenome.sh   # extract APOE seqs<br>bash 01.Variant_calling_PanGenome.sh<br># (then score haplotypes, see utils/)<br>Rscript 1.Fig3A.R … 6.Fig3F.R<br>``` | Diversity & ancestry panels |
| **4** | `Figure4/` | ```bash<br># ADNI cohort<br>cd Figure4<br>Rscript 0.EVO2_scoring_ADNI.R  # per-subject haplotypes<br>Rscript 1.Fig4BC.R 2.Fig4D.R 3.Fig4EF.R<br>``` | Diagnosis, cognition, PET |
| **5** | `Figure5/` | ```bash<br># UK Biobank scale<br>cd Figure5<br>Rscript 00.EVO2_scoring_UKB_VCF_to_haplotype.R<br>Rscript 01.EVO2_scoring_UKB_haplotype_to_sequence.R<br># MRI<br>Rscript 2.Fig5A_FAST.R 3.Fig5A_FAST_plot.R<br># SWI<br>Rscript 4.Fig5A_SWI.R 5.Fig5A_SWI_plot.R<br># Proteomics<br>Rscript 6.Fig5B_Olink.R 7.Fig5BC_Olink_plot.R<br># PheWAS<br>Rscript 8.Fig5D_PheWAS.R 9.Fig5D_PheWAS_plot.R<br>``` | MRI maps, volcano + enrichment, PheWAS Manhattan |

\*Adjust paths to your **Evo2 model**, **reference FASTA**, and dataset locations via environment variables (`APOE_DATA_DIR`, `EVO2_MODEL_DIR`) or script headers.

---

## 🗂️ Data & Model

| Resource | What it’s for | Access / Link |
| -------- | ------------- | -------------- |
| **Evo2-7B weights** | Sequence Δ-score generation | <https://github.com/oxford-deeplearning/evo2> |
| **GRCh38 ref (Chr 19 slice)** | Base sequence for variant/haplotype editing | Included in `utils/` or any GRCh38 FASTA |
| **ADNI WGS + phenotypes** | Individual haplotypes & PET / cognition (Fig 4) | Apply via <https://adni.loni.usc.edu> |
| **Human Pangenome assemblies (HPRC)** | Diverse haplotypes (Fig 3) | Public: <https://humanpangenome.org> |
| **UK Biobank genotypes, MRI, SWI, Olink, Phecodes** | Large-scale associations (Fig 5) | UKB Application; comply with MTA |
| **LDlinkR API (optional)** | LD calculations for plots | <https://ldlink.nci.nih.gov> |
| **R packages** | Stats & plots (`dplyr`, `ggplot2`, `robustbase`, `gprofiler2`, `cerebroViz`, …) | Auto-installed via `renv::restore()` |

> **Note** Some datasets (ADNI, UKB) require prior access approval.  
> Ensure credentials and local paths are configured before running scripts.


⸻

License

Pending – will follow journal policy (MIT for code, CC-BY 4.0 for docs).

⸻

Contact

Name	Email
Rui Zhu	rui.zhu.rz399@yale.edu
Xiaopu Zhou	xiaopu.zhou@sickkids.ca

**How this was shortened** — content reduced to ~50 %, with collapsible sections, tables, and direct code blocks so you can paste straight into `README.md` and get a clean, navigable document.
