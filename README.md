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

Reproduce the Figures

<details>
<summary><strong>▶ Figure 2 (Δ-score vs GWAS)</strong></summary>


cd Figure2
Rscript 0.EVO2_scoring_GWAS.R        # make SNP FASTA, call Evo2
Rscript 1.Fig2A.R                    # regional Manhattan plot
Rscript 2.Fig2B.R && 3.Fig2C.R       # Δ-score vs Z & β

Requires AD GWAS summary statistics (.tsv with chr, pos, p, beta).

</details>


<details>
<summary><strong>▶ Figure 3 (Pangenome)</strong></summary>


bash  Figure3/00.EVO2_scoring_PanGenome.sh   # extract APOE seq from assemblies
bash  Figure3/01.Variant_calling_PanGenome.sh
# score haplotypes (not included – see utils/evo2_score_batch.py)
Rscript Figure3/1.Fig3A.R  # ancestry-level plots
...

Needs minimap2, samtools, bcftools and HPRC assemblies.

</details>


<details>
<summary><strong>▶ Figure 4 (ADNI)</strong></summary>


cd Figure4
Rscript 0.EVO2_scoring_ADNI.R    # reconstruct & score haplotypes
Rscript 1.Fig4BC.R              # cohort overview
Rscript 2.Fig4D.R               # score vs cognition / CSF
Rscript 3.Fig4EF.R              # score vs amyloid PET

Provide ADNI WGS VCFs and phenotype CSVs.

</details>


<details>
<summary><strong>▶ Figure 5 (UK Biobank)</strong></summary>


cd Figure5
Rscript 00.EVO2_scoring_UKB_VCF_to_haplotype.R
Rscript 01.EVO2_scoring_UKB_haplotype_to_sequence.R
# MRI
Rscript 2.Fig5A_FAST.R  && 3.Fig5A_FAST_plot.R
# Proteomics
Rscript 6.Fig5B_Olink.R && 7.Fig5BC_Olink_plot.R
# PheWAS
Rscript 8.Fig5D_PheWAS.R && 9.Fig5D_PheWAS_plot.R

Requires access-approved UKB datasets (genotypes, imaging, Olink, phenotypes).

</details>



⸻

Data & Model

Resource	Access
Evo2-7B weights	https://github.com/oxford-deeplearning/evo2 (MIT license)
ADNI WGS & phenotypes	Application via ADNI data portal
Human Pangenome (HPRC)	Public assemblies: https://humanpangenome.org
UK Biobank	Application # ; abide by MTA
Reference genome	GRCh38 FASTA (chr19 slice shipped in utils/)


⸻

License

Pending – will follow journal policy (MIT for code, CC-BY 4.0 for docs).

⸻

Contact

Name	Email
Rui Zhu	rui.zhu.rz399@yale.edu
Xiaopu Zhou	xiaopu.zhou@sickkids.ca

**How this was shortened** — content reduced to ~50 %, with collapsible sections, tables, and direct code blocks so you can paste straight into `README.md` and get a clean, navigable document.
