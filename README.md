# Harnessing-DNA-Foundation-Models-for-Human-Population-Genomics-Insights-from-the-APOE-Locus
Harnessing DNA Foundation Models for Human Population Genomics – Code Repository

Overview

This repository contains the code accompanying the study “Harnessing DNA Foundation Models for Human Population Genomics: Insights from the APOE Locus.” In this work, we introduce a framework that integrates a DNA foundation model (Evo2) with population-scale genomic data to identify and prioritize functional genetic variation. We focus on the well-known APOE locus (chr19) – the strongest genetic risk region for Alzheimer’s disease – as a proof-of-concept to demonstrate how sequence-based deep learning models can complement traditional population genomics analyses.

Key aspects of the study and this code include:
	•	Evo2 model scoring: We utilize the Evo2-7B DNA foundation model (a large transformer trained on genome sequences) to score genetic variants and haplotypes in the APOE region. The model provides a Δ-score (delta score) for sequence variants/haplotypes, reflecting their predicted functional impact based on sequence context.
	•	Integration with genomic data: We apply model-based scores to real genomic datasets:
	•	Alzheimer’s Disease GWAS (ADNI/IGAP): to see if high model Δ-scores align with known AD association signals (Figure 2).
	•	Human Pangenome (HPRC) assemblies: to evaluate APOE haplotypes across diverse ancestries and detect functionally important differences (Figure 3).
	•	ADNI cohort whole genomes: to score individual APOE haplotypes and correlate with endophenotypes like cognitive tests and amyloid PET imaging (Figure 4).
	•	UK Biobank: to score APOE haplotypes in a biobank-scale population and perform phenome-wide association analyses, including brain imaging, proteomics, and clinical phenotypes (Figure 5).
	•	Haplotype-level analysis: Because the APOE region contains extended haplotypes (e.g., the ε4-bearing haplotype spans ~80 kb across APOE, TOMM40, APOC1, etc.), our approach reconstructs full haplotype sequences for individuals and scores them holistically with the model. This enables per-haplotype risk profiling beyond single-variant tests.

Overall, this repository provides all scripts needed to reproduce Figures 2–5 of the paper, which illustrate how model-derived sequence scores can highlight likely functional variants and haplotypes and how these scores correlate with population genetic signals and disease phenotypes.

Repository Structure

The code is organized by figure, with separate folders for Figure 2, 3, 4, and 5 analyses. Each folder contains scripts (R or shell scripts) used to generate the results and visualizations for that figure. Below is an overview of the repository contents:

.
├── Figure2
│   ├── 0.EVO2_scoring_GWAS.R
│   ├── 1.Fig2A.R
│   ├── 2.Fig2B.R
│   ├── 3.Fig2C.R
│   ├── APOE_ref_demo.txt
│   ├── SNP_demo.txt
│   └── results
│   
├── Figure3
│   ├── 00.EVO2_scoring_PanGenome.sh
│   ├── 01.Variant_calling_PanGenome.sh
│   ├── 1.Fig3A.R
│   ├── 2.Fig3B.R
│   ├── 3.Fig3C.R
│   ├── 4.Fig3D.R
│   ├── 5.Fig3E.R
│   └── 6.Fig3F.R
│   
├── Figure4
│   ├── 0.EVO2_scoring_ADNI.R
│   ├── 1.Fig4BC.R
│   ├── 2.Fig4D.R
│   └── 3.Fig4EF.R
│   
└── Figure5
    ├── 00.EVO2_scoring_UKB_VCF_to_haplotype.R
    ├── 01.EVO2_scoring_UKB_haplotype_to_sequence.R
    ├── 1.TableS10.R
    ├── 2.Fig5A_FAST.R
    ├── 3.Fig5A_FAST_plot.R
    ├── 4.Fig5A_SWI.R
    ├── 5.Fig5A_SWI_plot.R
    ├── 6.Fig5B_Olink.R
    ├── 7.Fig5BC_Olink_plot.R
    ├── 8.Fig5D_PheWAS.R
    └── 9.Fig5D_PheWAS_plot.R

	•	Figure2/ – Code for analyzing GWAS signals vs. Evo2 model scores in the APOE region (Alzheimer’s disease GWAS integration).
	•	Figure3/ – Code for analyzing APOE region haplotypes from the Human Pangenome Project and diverse ancestries.
	•	Figure4/ – Code for analyzing individual APOE haplotypes in the ADNI cohort and correlating model scores with AD-related phenotypes.
	•	Figure5/ – Code for large-scale analysis in UK Biobank: haplotype identification, model scoring, and associations with imaging, proteomic, and clinical phenotypes.

Each script is prefixed by a number indicating its role or order in the analysis. For example, 0.* or 00.* scripts typically perform data processing or scoring (often generating intermediate data like model scores), while subsequent numbered scripts (1, 2, 3, …) perform plotting or statistical analysis for specific figure panels.

Requirements

To run these analyses, you will need the following software and data resources:
	•	R environment: The R scripts were tested on R 4.x and require various R packages. Key packages include:
	•	Data manipulation and visualization: dplyr, data.table, tidyr, ggplot2, ggrepel, ggridges, ggsci, colorspace.
	•	Statistical analysis: MethComp (for Deming regression), brms and posterior (for Bayesian modeling), ppcor (partial correlations), robustbase (robust regression), lme4 (mixed models, if needed), RNOmni (rank normalization).
	•	Genomic/association analysis: LDlinkR (for LD calculations or regional plots), topr (for LocusZoom-style plots), vcfR (VCF file parsing).
	•	Enrichment and annotation: gprofiler2 (GO/KEGG enrichment analysis).
	•	Domain-specific: cerebroViz (for brain region visualizations on MRI).
	•	(See the individual library(...) calls at the top of each script for the full list of packages used. Ensure these are installed before running.)
	•	Command-line tools: Some analyses (especially in Figure3) require external genomic tools:
	•	minimap2 (for sequence alignment of assemblies to reference).
	•	samtools and bcftools (for processing alignments and variant calling).
	•	These should be installed and available in your $PATH when running the shell scripts.
	•	Data resources:
	•	Reference genome sequence: A FASTA sequence of the APOE locus region (human chr19) is needed for aligning and variant calling. In our pipeline we used GRCh38 coordinates approximately spanning 44,841,743–44,921,743 on chr19 (which covers APOE and adjacent genes). An example reference file APOE_ref_demo.txt is provided as a template (for demonstration) in Figure2, but you should replace it with the actual reference sequence for full analysis.
	•	GWAS summary statistics: The AD GWAS results (e.g., from IGAP or other studies) for the APOE region are expected for plotting Figure 2A and computing correlations in 2B/2C. The code expects summary stats including SNP IDs, p-values, effect sizes, etc.
	•	Human Pangenome Project assemblies: Figure3’s pipeline will download the human pangenome graph assembly (GFA) for chromosome 19 from the Human Pangenome Reference Consortium. Ensure you have a fast internet connection and sufficient storage for these assemblies. The script will automatically download needed files (if not already present) to a specified path.
	•	ADNI WGS data: Figure4 analysis requires individual-level variant data in the APOE region for participants of the Alzheimer’s Disease Neuroimaging Initiative (ADNI). This is typically provided as VCFs (the code expects VCF files per sample or a merged VCF) and associated phenotype data:
	•	Genotypes: The scripts parse VCF to reconstruct each individual’s two haplotype sequences.
	•	Phenotypes: ADNI metadata such as diagnosis, age, sex, APOE genotype, cognitive scores, and imaging (amyloid PET SUVR) should be available. The code will attempt to merge haplotype scores with these data for regression analyses.
	•	UK Biobank data: Figure5 is based on UK Biobank participants.
	•	Genotypes: Access to UKB imputed or WGS data is required to extract APOE region haplotypes. The pipeline (00 and 01 scripts) expects a phased VCF or similar haplotype data for the APOE region of UKB participants.
	•	Phenotypes: Various UKB phenotypes are used:
	•	Brain MRI measures: Regional gray matter volumes (from structural MRI) and SWI (Susceptibility-Weighted Imaging) metrics. These need pre-processing (e.g., volumetric measures normalized for intracranial volume, SWI measures of iron or microbleeds). The code assumes these measures are available for each individual in a tabular format.
	•	Proteomics (Olink) data: If analyzing protein associations (Fig5B/C), Olink proteomic assay data for UKB participants (protein levels) should be available.
	•	Clinical diagnoses/traits: For the PheWAS (Fig5D), a broad set of phenotypes or disease outcomes (possibly encoded as PheCodes or ICD-derived traits) are used. Ensure you have a prepared phenotype table with case-control status for many traits, along with a mapping to categories (for coloring the Manhattan plot).
	•	Covariates: Typically, analyses will adjust for covariates like age, sex, principal components of ancestry, etc. The code might expect these in the phenotype data merges.
	•	Local ancestry data (optional): In Figure3C, local ancestry inference results from RFMix2 were used to correlate ancestry proportions with model scores along the region. If you wish to reproduce that analysis, you should have run RFMix2 on APOE region for relevant samples (e.g., in a diverse dataset or the 1000 Genomes samples) to produce local ancestry calls. The script will read and integrate those results if available.

Note on Evo2 model: The Evo2 DNA foundation model (7B parameters) is a core part of this analysis. You will need access to the model to score sequences:
	•	In our study, Evo2 was used in a zero-shot manner to obtain a “functional score” for sequences. Practically, this means for each sequence (reference or mutated haplotype), the model assigns a likelihood or a functional prediction, and we derive a Δ-score representing the difference in log-likelihood or functional score caused by a variant/haplotype.
	•	Ensure you have the capability to run the Evo2 model (either via an API or locally with the model weights) since the scripts will either call out to the model or expect precomputed scores. (The specifics of model scoring may require additional code or configuration not included in this repository — for example, loading the model in Python and querying sequences. You may need to integrate that step manually or use provided placeholders where applicable.)

Finally, before running any scripts, adjust file paths and input filenames as needed. Many scripts have configuration sections at the top (for example, the shell scripts in Figure3 define paths to tools and data directories, and R scripts may assume working directory structure as given in this repo). Update those paths according to your environment.

Figure 2 – Variant Scoring and GWAS Integration

Figure 2 in the paper examines how the Evo2 model’s sequence-based scores correlate with Alzheimer’s GWAS signals in the APOE region. The code in Figure2/ accomplishes the following main tasks:
	•	Preparing variant sequences and model scores: We take a list of SNPs in the APOE region and generate mutant sequences for each SNP (substituting the alternate allele into the reference sequence). Evo2 is used to score these sequences, and the difference from reference (Δ-score) is computed.
	•	Comparing to GWAS results: We then compare the model’s Δ-scores for variants to GWAS statistics (–log10 p-values, effect sizes, Z-scores) to see if high-impact model predictions coincide with strong association signals.
	•	Visualization: The scripts produce LocusZoom-style regional plots and scatter plots showing the relationship between model scores and GWAS metrics.

Scripts for Figure 2:

Script	Purpose
0.EVO2_scoring_GWAS.R	Variant processing and sequence scoring. Loads the APOE reference sequence (example provided in APOE_ref_demo.txt) and a list of SNPs (SNP_demo.txt). It verifies reference alleles, filters out any mismatches, and then creates mutated FASTA sequences for each SNP’s alternate allele. These sequences are scored with the Evo2 model to compute Δ-scores. Outputs: filtered SNP list with correct ref/alt alleles, matched SNP list (reference vs alternate allele sequences), and a set of mutated sequence files (in results/FASTA_modified/) with their model scores.
1.Fig2A.R	Regional association plots. Generates two LocusZoom-style plots for the APOE region: (1) an AD GWAS Manhattan plot (–log₁₀p vs genomic position) highlighting the lead SNP rs429358 (APOE-ε4), and (2) an Evo2 Δ-score plot (model-predicted impact vs position) highlighting rs7412 (APOE-ε2). These plots illustrate the correspondence (or differences) between statistical association peaks and model-predicted functional variants.
2.Fig2B.R	Score vs. GWAS scatter plots (all variants). Plots the relationship between Evo2 Δ-scores and GWAS effect estimates for all SNPs in the region. It creates two scatter plots: (1) Δ-score vs. GWAS Z-score and (2) Δ-score vs. GWAS β (effect size). Each plot includes a fitted linear regression line to assess correlation. This helps evaluate whether variants the model predicts as highly functional also show large effect sizes or strong significance in the GWAS.
3.Fig2C.R	Score vs. GWAS scatter plots (clumped lead SNPs). Similar to Fig2B, but restricted to a set of 28 independent lead SNPs (obtained via clumping or fine-mapping of the region). By focusing on the top signals, these two scatter plots (Δ-score vs Z, and Δ-score vs β) test if the model better predicts the most significant causal variants. A stronger correlation here would suggest the model is effectively prioritizing true causal hits.

Usage notes: To reproduce Figure 2 results, run the scripts in order. First execute 0.EVO2_scoring_GWAS.R to generate the necessary data (it will produce the mutated sequences and a table of scores). Then run 1.Fig2A.R to create the regional plots, followed by 2.Fig2B.R and 3.Fig2C.R for the scatterplots. Ensure you have the GWAS summary stats file ready (with at least SNP positions, p-values, etc.) as the plotting scripts will likely read those results. The output plots can be saved as PDF or image files as configured in the scripts (look for ggsave() or similar commands at the end of each script to specify output filenames if needed).

Figure 3 – Human Pangenome Haplotype Analysis

Figure 3 explores APOE region haplotypes across diverse human populations using the Human Pangenome Reference Consortium (HPRC) data and population genetics approaches. The goal is to see how the Evo2 model scores vary among haplotypes from different ancestries and whether those patterns reflect known differences (e.g., the attenuated effect of APOE-ε4 in African ancestry).

This analysis involves a multi-step pipeline:
	1.	Retrieving APOE sequences from the Pangenome: We leverage assemblies from the HPRC (which includes genomes across global populations). We extract the APOE locus sequences from each assembly for comparison.
	2.	Variant calling on assemblies: Using the extracted sequences, we identify sequence variants in the APOE region (relative to the reference genome). This gives a comprehensive set of polymorphisms present across diverse lineages.
	3.	Scoring haplotypes: The unique haplotype sequences from these assemblies (and potentially from 1000 Genomes or other reference panels) are scored by the Evo2 model to obtain Δ-scores for each haplotype.
	4.	Population analyses: Using the scored haplotypes, we perform several analyses:
	•	Correlate model scores with APOE gene expression differences and known effect sizes by ancestry.
	•	Examine haplotype-specific correlations, including distinguishing the classic ε2/ε3/ε4 alleles.
	•	Analyze local ancestry patterns: how segments of the APOE region inherited from different ancestral backgrounds might influence the model’s score.
	•	Compare model predictions with observed risk (odds ratios) of APOE-ε4 in different populations.

Scripts for Figure 3:

Data Processing (Shell pipeline):

Script	Purpose and steps
00.EVO2_scoring_PanGenome.sh	Pangenome assembly processing pipeline. This is a Bash script that automates obtaining and processing the human pangenome assemblies for the APOE region. It performs the following steps:  • Download the chromosome 19 GFA assembly (graph-based assembly) from the HPRC public repository.  • Align the APOE reference sequence to the pangenome using minimap2, producing PAF alignment files for each assembly.  • Concatenate and parse alignment files to identify aligned segments.  • Extract the APOE locus sequence from the reference genome (for coordinates ~44.84–44.92 Mb on chr19) and the corresponding sequences from each assembly’s alignment.  Output: a set of APOE region sequences for each assembly/genome in FASTA format (or a similar consolidated file), to be scored by Evo2. (Paths and parameters within should be edited to match your system; e.g., location of minimap2 and output directories.)
01.Variant_calling_PanGenome.sh	Variant calling on assemblies. This script takes the extracted APOE region sequences from all assemblies and identifies genetic variants:  • Uses samtools and bcftools to create a multi-sequence alignment (via SAM/BAM) and call variants relative to the reference.  • Outputs a VCF file containing all SNPs and indels in the APOE region across the pangenome assemblies.  These variants represent population-diverse variation in APOE and will be used for downstream analysis and haplotype classification. (Dependencies: samtools, bcftools; ensure paths in the script are correct.)

After running the above, you should have: (a) a collection of APOE sequences from diverse genomes, and (b) a VCF of variants. Next steps (not fully scripted here) would include scoring each unique haplotype sequence with Evo2. In practice, one would take each assembly’s APOE sequence and compute its Δ-score using the model (comparing to the reference sequence). The resulting scores can then be associated with haplotype features (like whether it carries ε4, etc.). The R scripts below assume such scoring has been done and relevant summary data prepared (e.g., mean scores by population, etc.).

Analysis and Plotting (R scripts):

Script	Purpose and analysis details
1.Fig3A.R	Ancestry-stratified APOE effect vs. model score. Calculates and visualizes the relationship between Evo2 scores and APOE functional effects across major ancestral groups. It takes the mean Δ-score for APOE haplotypes in each super-population (AFR, AMR, EAS, EUR, SAS) and compares it to a measure of APOE’s effect (for example, an expression level or effect size of APOE variants in that population). A Bayesian Deming regression (accounting for errors in both variables) is fit to these points. The script outputs the correlation and Bayesian p-values, and generates a scatter plot with a regression line, where each point is an ancestry group (color-coded) with error bars. This addresses whether ancestry differences in APOE’s impact (such as lower ε4 effect in Africans) align with model-predicted functional differences.
2.Fig3B.R	Per-haplotype expression vs. score correlation. Merges haplotype-specific APOE expression data with their model Δ-scores. Each data point here is a specific APOE haplotype (for instance, distinct haplotypes in the pangenome or 1000G data). Haplotype are annotated as “leading” vs “alternative” (which could refer to the major vs minor haplotype for a given individual or population). The script then plots the correlation between a haplotype’s Evo2 score and its APOE gene expression level (or effect on expression), highlighting whether the primary haplotype carries the majority of expression. This yields a scatter plot assessing if the model’s functional score correlates with actual gene expression differences across haplotypes.
3.Fig3C.R	Local ancestry vs. functional score (partial correlations). Integrates RFMix2 local ancestry output with model scores along the APOE region. The APOE locus may be divided into sub-regions or SNP bins, and for each region the script examines the relationship between the fraction of a particular ancestry and the Evo2 Δ-score, while controlling for other ancestries (partial Spearman correlation). Essentially, for each ancestry X and region segment Y, it computes how the presence of ancestry X in that segment correlates with the haplotype’s score, controlling for global effects. The result is a set of correlation coefficients (and significance) for each ancestry–region combination, which can be output as a heatmap or table. This helps identify if certain ancestry segments (e.g., an African-derived segment within an otherwise European haplotype) contribute to lower or higher predicted functional impact.
4.Fig3D.R	Haplotype classification by ε2/ε3/ε4 and score distribution. Classifies each haplotype based on the presence of APOE ε4 or ε2 alleles (using SNPs rs429358 and rs7412 status). For example, haplotypes are labeled as ε4-carrying, ε2-carrying, or neutral (ε3/ε3 background). The script then compares the Evo2 Δ-scores across these categories – e.g., do ε4 haplotypes have systematically higher (more damaging) scores than ε3? It likely produces a plot (such as a boxplot or density plot) of the distribution of model scores for haplotypes in each APOE allele category. This illustrates whether the model recognizes the higher risk of ε4 at the sequence level.
5.Fig3E.R	ε4 haplotypes – score distribution by ancestry. Focuses on haplotypes that carry the ε4 allele. It merges the model Δ-scores for those haplotypes with sample metadata (particularly ancestry information), and then visualizes the distribution of scores for ε4-haplotypes across different ancestral groups. For instance, it may show whether ε4 haplotypes found in African-ancestry individuals have different model-predicted impact than ε4 haplotypes in European or Asian individuals. This analysis relates to observations that the ε4 risk is modulated by genetic background – the model might score an African ε4 haplotype as less damaging if surrounding sequence context is protective. The output could be a ridge plot or boxplot stratified by ancestry, for ε4 carriers only.
6.Fig3F.R	Model score vs. observed ε4 risk across populations. Compares the normalized Evo2 Δ-score for ε4 haplotypes in each population to the actual odds ratio (OR) of APOE ε4 for Alzheimer’s disease in that population. For example, if European ε4 has OR ~3, African ε4 OR ~1.5, etc., and the model scores for those haplotypes might differ, we check if there’s a correlation. The script merges a table of ε4 ORs by ancestry (from epidemiological studies) with the haplotype score data, and plots them with error bars (for OR confidence intervals, etc.). It produces a correlation plot where each point is an ancestry group, x-axis being model score and y-axis the observed OR (or vice versa), showing how well the model’s predictions of functional impact align with real-world risk.

Usage notes: The shell scripts 00.EVO2_scoring_PanGenome.sh and 01.Variant_calling_PanGenome.sh should be run first to generate the APOE sequences and variant calls from the pangenome data. Adjust paths for tools and output directories as needed. After running them, ensure you have performed Evo2 scoring on the assembled sequences (this step may not be fully automated in the script – you might need to run a separate Python script or command to feed the sequences to the Evo2 model and collect scores). Once you have the model Δ-scores and any other necessary data (like expression levels or ORs by population), you can proceed with the R scripts 1.Fig3A.R through 6.Fig3F.R. These may require assembling input tables from the previous steps (for instance, a table of haplotypes with their scores, ancestry assignments, and any phenotypes). Check each script’s top section for expected input file names or data formats:
	•	For example, 1.Fig3A.R expects a small table of mean scores and sem (standard errors) by superpopulation.
	•	3.Fig3C.R expects RFMix local ancestry output (perhaps as a matrix of ancestry proportions per region).
	•	6.Fig3F.R expects a file or hardcoded values for ε4 odds ratios by ancestry.

Make sure to provide these inputs in the correct format. Running each script will produce the panels for Figure 3 (scatter plots, correlation heatmaps, boxplots, etc., as described above). Save the plots to files as needed (the scripts likely contain code to save plots; modify the paths if you want them in a specific directory).

Figure 4 – ADNI Cohort Haplotype Analysis

Figure 4 shifts focus to individual-level analysis in a well-phenotyped cohort: the Alzheimer’s Disease Neuroimaging Initiative (ADNI). Here, the goal is to assess whether the Evo2 model’s haplotype scores correlate with clinical and biomarker data in a real patient cohort. We use whole-genome sequences from ADNI to extract each participant’s APOE haplotypes, score them, and then examine associations with Alzheimer’s-related phenotypes (diagnosis, cognitive scores, imaging biomarkers).

Summary of Figure 4 analysis:
	•	Haplotype reconstruction and scoring: For each ADNI subject, reconstruct the two APOE region haplotypes from their WGS VCF, substitute any variants into the reference sequence to get the exact sequence, and obtain a model Δ-score for each haplotype.
	•	Clinical data merge: Combine the haplotype score data with ADNI participants’ metadata (such as APOE genotype, diagnosis status, age, etc.) and downstream phenotypes.
	•	Associations with phenotypes: Analyze whether a participant’s haplotype score (particularly the more pathogenic haplotype’s score) is associated with:
	•	Diagnosis of AD/MCI vs control (baseline diagnostic status, or conversion status).
	•	Endophenotypes: e.g., cognitive test scores, CSF biomarkers, etc. (In our code, we specifically look at baseline cognitive or biomarker measurements.)
	•	Imaging biomarkers: here, an emphasis on amyloid PET burden (measured by PET SUVR in regions like frontal cortex, etc. – a hallmark AD biomarker).
	•	Visualization: Produce figures that show differences or correlations, such as bar charts for diagnosis groups, scatterplots of haplotype score vs. phenotype measure, including regression lines or confidence intervals.

Scripts for Figure 4:

Script	Purpose and analysis steps
0.EVO2_scoring_ADNI.R	APOE haplotype reconstruction and scoring for ADNI. This script reads in the APOE reference sequence (similar to previous, a FASTA-like text for the region) and parses ADNI whole-genome VCF data for the APOE locus. It iterates through individuals, and for each individual:  • Determines their genotype at each variant in the region, splitting into two haplotypes (phasing info from VCF if available, or by reference phasing if needed).  • Reconstructs the haplotype DNA sequences by introducing the individual’s variants into the reference sequence for each chromosome copy.  • Scores each haplotype sequence using the Evo2 model (likely via an internal call or an external script – ensure the model is accessible).  • Records the results in a table, including subject ID, haplotype identifier (1 or 2), and the model Δ-score.  Output: A dataset of ADNI haplotypes with their scores. This is the foundation for subsequent analysis scripts.
1.Fig4BC.R	Demographics and genotype summaries (Fig 4B,4C). Merges the haplotype score data with ADNI diagnosis and demographic metadata. It likely classifies individuals by their APOE genotype (ε2/ε3/ε4 status) and perhaps by whether their haplotype score is high or low. The script generates summary statistics and plots for figure panels 4B and 4C. For instance:  • Figure 4B: could be a bar chart or table summarizing the ADNI cohort by APOE genotype and diagnosis (showing sample counts, ages, etc., possibly to illustrate the makeup of the cohort).  • Figure 4C: possibly a visualization of how the haplotype Δ-scores distribute among different diagnostic groups or APOE genotypes (e.g., comparing scores of ε4 carriers vs non-carriers, or AD vs control). This script ensures a baseline understanding of the data and that haplotype scoring aligns with known factors (e.g., ε4 carriers should have higher scores on one haplotype).
2.Fig4D.R	Haplotype score vs baseline endophenotypes (Fig 4D). Integrates the haplotype scores with baseline clinical/endophenotype data from ADNI. This may include cognitive test scores (e.g., MMSE, ADAS-Cog) or fluid biomarkers (like CSF Aβ42, tau) at baseline. The script then performs association analysis – likely a robust linear regression (using lmrob from robustbase) – to test if an individual’s haplotype Δ-score is a significant predictor of these baseline measures. Covariates (age, sex, education, etc.) might be included. The result for Figure 4D is likely a plot or set of results showing a significant correlation; for example, individuals with higher (more damaging) haplotype scores might have lower cognitive performance or abnormal biomarker levels. The script could output a regression summary and create a scatter plot with a regression line illustrating the relationship for a particular endophenotype that showed significance.
3.Fig4EF.R	Haplotype score vs amyloid PET (Fig 4E, 4F). Focuses on imaging data, specifically amyloid PET SUVr (a measure of amyloid plaque deposition in the brain). The script merges the haplotype score data with PET measurements at baseline (for regions of interest or a composite measure). It then runs linear models to see if haplotype Δ-score predicts amyloid burden. Figure 4E and 4F likely present this analysis in two ways:  • Fig 4E: possibly a scatter plot of haplotype score vs PET SUVr (with a fit line), showing whether higher scores correlate with higher amyloid.  • Fig 4F: perhaps a grouped comparison, e.g., dividing subjects into tertiles of score or by APOE genotype and comparing mean amyloid levels, illustrating how those with high model-predicted risk haplotypes accumulate amyloid.  The script outputs the statistical results (p-values, coefficients) and generates the plots for these panels. If robust regression or adjustments are used (to downweight outliers or adjust for covariates like age), those details are handled here.

Usage notes: Start with 0.EVO2_scoring_ADNI.R to produce the core haplotype score data for all ADNI subjects. This requires having all individual VCF files (or a combined VCF) for the APOE locus and may be computationally intensive if done for many individuals. After obtaining the haplotype score table, ensure you have ADNI phenotype data files (such as diagnosis info, cognitive scores, biomarker levels, PET SUVR data) readily accessible – these should be loaded or read by the subsequent scripts. You may need to modify file names/paths in scripts 1–3 to point to your phenotype data files. Once the data is in place:
	•	Run 1.Fig4BC.R: it will likely output some summary and save plots for B and C. Check the script for any file outputs (it might print a table or save plots).
	•	Run 2.Fig4D.R: this will perform regressions; ensure any specific phenotype names or columns used in the script match your data. It should produce output for panel D.
	•	Run 3.Fig4EF.R: ensure the PET SUVR data is correctly linked. This will yield results for panels E and F.

Each of these scripts may save plots to files (look for ggsave or similar in the code). If needed, adjust the output filenames or paths. Interpret the results in context: a significant finding in these analyses would support the hypothesis that the model’s haplotype scores reflect biologically meaningful differences (e.g., a high-scoring haplotype corresponds to higher amyloid in the brain).

Figure 5 – UK Biobank Haplotype Associations

Figure 5 extends the analysis to the UK Biobank (UKB), a much larger dataset, to demonstrate scalability and to perform a broad phenome-wide evaluation of APOE haplotype scores. The idea is to use the Evo2 model to score all distinct haplotypes in the APOE region found in ~500,000 UKB participants, then test how those scores associate with a variety of phenotypes, including neuroimaging measures, proteomics, and clinical diagnoses.

Key components of Figure 5 analysis:
	•	Haplotype identification in UKB: Because of the large sample size, rather than scoring every individual’s two haplotypes separately, the approach is likely to identify unique haplotypes present in the population. UKB data (imputed or WGS) can be processed to list unique haplotype sequences or at least unique combinations of key variants (like tagging SNPs defining haplotypes).
	•	Haplotype sequence reconstruction and scoring: For each unique haplotype (or each individual if needed), construct the sequence of the APOE region and obtain its Evo2 Δ-score. Given UKB’s scale, this likely involved focusing on the most common haplotypes and those carrying interesting variants (like ε4, ε2, and other combinations).
	•	Phenotype associations: Using the scores (especially assigning each individual a haplotype score – possibly the “leading” haplotype score meaning the higher of their two haplotype scores), perform association tests with:
	•	Brain MRI traits: particularly structural MRI volumes (gray matter volume in various brain regions) and Susceptibility-Weighted Imaging (SWI) metrics (which relate to iron deposition or microbleeds, potentially relevant in AD).
	•	Plasma proteomics: levels of hundreds of proteins (Olink assay) to see if APOE-related pathways (like cholesterol, inflammation proteins) show association with haplotype score.
	•	Clinical outcomes (PheWAS): a wide array of diagnoses or quantitative traits to see if any phenotype beyond AD is associated with the model score (for example, other neurological diseases, cardiovascular traits, etc., reflecting pleiotropy of APOE region).
	•	Visualization: Fig 5 includes multiple panels:
	•	5A likely covers MRI results (possibly two sub-parts: one for structural MRI, one for SWI).
	•	5B and 5C cover proteomic results (a volcano plot of protein associations and an enrichment analysis of significant proteins).
	•	5D is a PheWAS Manhattan plot summarizing associations across many phenotypes, highlighting where APOE haplotype score shows significant effects (AD should appear here, validating the approach).

Scripts for Figure 5:

Data Processing:

Script	Purpose and process
00.EVO2_scoring_UKB_VCF_to_haplotype.R	Parse UKB genotypes to haplotypes. This R script reads UK Biobank genetic data for the APOE region (likely a VCF or PLINK data for chr19 around APOE) and identifies haplotypes. Steps include:  • Read in genotype data for all individuals (could be large; possibly done in chunks).  • If phased data is available, directly extract haplotype sequences (if not, phasing might be required prior).  • Determine unique haplotypes and assign each individual’s haplotypes an ID or representation. This might involve focusing on key variant combinations (e.g., tagging SNPs to differentiate major haplotypes) to reduce complexity.  • Output: a mapping of individuals to haplotype IDs, and a list of unique haplotype allele compositions to be reconstructed.
01.EVO2_scoring_UKB_haplotype_to_sequence.R	Haplotype sequence assembly and scoring. Takes the unique haplotype definitions from the previous step and builds actual DNA sequences for the APOE region for each unique haplotype. This likely involves starting from the reference sequence and inserting the variants that define that haplotype. It then uses the Evo2 model to score each haplotype sequence (again, obtaining a Δ-score for each).  • The script may loop through haplotypes, create FASTA entries, and call the model (possibly via a system call to a Python script or an API).  • Output: a table of haplotypes with their Evo2 scores. Additionally, each individual in UKB can now be assigned a score_leading (the higher of their two haplotype scores) and perhaps a score_secondary (the lower). The “score_leading” is often used to represent an individual’s genetic risk (since, for example, an ε4/ε3 genotype individual’s risk is driven by the ε4 haplotype primarily).

After these two steps, you should have: (a) a list of unique haplotypes with scores, and (b) each UKB participant annotated with the score of their highest-risk haplotype (and possibly other haplotype info). The remaining scripts perform analyses on this data.

Analysis and Visualization:

Script	Purpose and details
1.TableS10.R	Logistic regression for AD (Supplementary Table 10). This script integrates the UKB haplotype score data with phenotype data to specifically examine Alzheimer’s disease status. It likely pulls the AD case-control status for UKB participants (e.g., ICD codes or self-reported diagnosis of AD or proxy dementia diagnoses) and performs a logistic regression with Evo2 score as a predictor, adjusting for covariates (age, sex, principal components, APOE genotype perhaps). The output would be an odds ratio for haplotype score predicting AD risk, along with confidence interval and p-value, which is reported as Table S10 in the paper. This serves as a direct validation that the model score is associated with actual AD cases in the biobank.
2.Fig5A_FAST.R	Association with structural MRI (FAST) – Analysis. This script tests the association between haplotype scores and brain gray matter volumes measured by MRI. “FAST” here refers to a set of brain regions or an analysis pipeline (it might be an acronym or just a label we use for the structural MRI analysis). The script:  • Merges each individual’s score (likely the leading haplotype score) with their MRI phenotypes (volumes of various brain structures, normalized for head size).  • Performs robust linear regression for each brain region’s volume ~ score (controlling for covariates like age, sex, MRI site, etc., as appropriate).  • Outputs the regression coefficients and p-values for the association between haplotype score and volume in each region. It also applies multiple testing correction (FDR) across regions.  These results indicate if higher genetic risk scores correlate with structural brain atrophy in specific areas.
3.Fig5A_FAST_plot.R	Structural MRI results – Plotting. Takes the results from 2.Fig5A_FAST and generates visualizations. Likely it uses the cerebroViz package to create brain maps highlighting regions where associations were found. Specifically:  • Annotates brain regions with their association statistics (e.g., which regions had significant volume reduction per unit increase in score).  • Creates lateral and medial brain hemisphere plots indicating the significant regions (perhaps coloring regions by the regression coefficient or –log₁₀ p-value).  • Possibly separates results by left/right hemisphere or region type.  The output is Figure 5A (first part), showing a brain diagram with affected regions, illustrating a structural phenotype associated with the APOE haplotype score.
4.Fig5A_SWI.R	Association with SWI MRI metrics – Analysis. This script focuses on Susceptibility-Weighted Imaging (SWI) outcomes. SWI in UKB could include measures like iron deposition in deep brain nuclei or microhemorrhages count. The script:  • Merges haplotype scores with SWI-derived metrics for participants.  • Filters and normalizes these metrics as needed (some may need log-transform or rank normalization due to skewness).  • Performs robust regression or linear models for each SWI metric ~ score (with covariates).  • Outputs coefficients and p-values, determining if haplotype score predicts any SWI measure (for instance, higher score correlating with more iron in certain brain areas).
5.Fig5A_SWI_plot.R	SWI results – Plotting. Similar to the FAST plotting, this script visualizes the SWI association results. If the SWI measures correspond to specific brain regions or structures, it may again use brain plots or bar charts:  • Highlights any significant associations found (e.g., perhaps increased iron in basal ganglia with higher score, etc.).  • Could produce a brain map focused on those specific SWI metrics or a summary figure.  Combined with the FAST results, parts of Figure 5A will illustrate how the haplotype score relates to different brain imaging phenotypes (structural atrophy and microhemorrhage/iron deposition patterns).
6.Fig5B_Olink.R	Association with plasma proteomics – Analysis. This script tests the haplotype score against a broad panel of protein biomarkers measured by the Olink platform in UKB. It likely:  • Reads in normalized protein level data for a set of proteins (hundreds of them) for individuals that have both genetic data and proteomic data.  • For each protein, fits a regression or correlation model: protein level ~ haplotype score + covariates (age, sex, principal components, possibly technical covariates from the proteomic assay).  • Collects the beta coefficients and p-values for the association of score with each protein.  • Applies multiple testing correction to identify significant protein associations.  This yields a list of proteins whose levels are significantly associated with the APOE haplotype score. These might include, for example, proteins related to lipid metabolism or neurodegeneration. (This script may not itself produce a figure, but prepares the data for Fig5B/5C.)
7.Fig5BC_Olink_plot.R	Proteomics results – Volcano plot & enrichment (Fig 5B, 5C). Using the output from 6.Fig5B_Olink, this script creates:  • Figure 5B: a volcano plot of all proteins, where the x-axis is the effect size (or direction of association) of haplotype score on protein level, and the y-axis is –log₁₀ p-value. Each point is a protein, and significant hits (after FDR correction) are highlighted. This plot quickly shows if any proteins stand out as strongly associated (either increased or decreased with higher genetic risk score).  • Figure 5C: a pathway enrichment analysis of the significant proteins. The script uses gprofiler2 to perform GO and KEGG pathway enrichment on the set of proteins that were significant in the association test. The results are then plotted, perhaps as a bar chart of top enriched pathways or processes. This helps interpret the proteomic hits in a biological context (e.g., enriched pathways might include cholesterol transport or inflammation, aligning with APOE’s known functions).  The combination of B and C panels demonstrates not only specific protein biomarkers linked to the haplotype score, but also the broader biological themes of those biomarkers.
8.Fig5D_PheWAS.R	Phenome-wide association study (PheWAS) – Analysis. This script casts a wide net by testing the haplotype score against many phenotypes in UKB (beyond AD). It likely works with a pre-defined list of phenotype codes (PheCodes or UKB trait IDs) that cover various disease categories. Steps:  • For each phenotype (e.g., ischemic heart disease, stroke, diabetes, etc.), perform a logistic regression (for case-control outcomes) or linear regression (for quantitative traits) with haplotype score as a predictor (plus covariates).  • Record the effect size and p-value for each association.  • Adjust for multiple testing (e.g., use an FDR threshold or Bonferroni for significance).  • Organize results by phenotype category (like grouping by cardiovascular, neurological, metabolic, etc., if using PheCodes).  The output is a comprehensive list of traits and their association statistics with the model score. We expect Alzheimer’s disease or dementia to be among the top hits (validating the approach), but this analysis may reveal other associations (for instance, APOE is also known to affect lipid levels, so we might see signals in cholesterol-related traits, etc.).
9.Fig5D_PheWAS_plot.R	PheWAS results – Manhattan plot (Fig 5D). This script visualizes the results from 8.Fig5D_PheWAS as a Manhattan plot. Instead of genomic positions, the x-axis here is different phenotypes grouped by category (each category can be plotted in a separate section along the axis), and the y-axis is the –log₁₀ p-value of the association. The script:  • Assigns each phenotype a category and a color (e.g., all cardiovascular outcomes in one color, all neurological in another, etc.).  • Plots a point for each phenotype’s association p-value. Horizontal lines indicate significance thresholds (e.g., FDR cutoff).  • Labels or highlights those phenotypes that reach significance. We expect to see Alzheimer’s disease on the plot as a genome-wide significant hit (since the haplotype score was designed around AD), perhaps along with other significant phenotypes (maybe lipid disorders or longevity, etc., if any).  This figure thus gives a broad view of what traits are affected by APOE-region variation as captured by the model’s score.

Usage notes: The UKB analysis is the most data-intensive. Ensure that you have sufficient memory and that any large data (genotypes, imaging, proteomics) are accessible in your R environment. A suggested order of execution:
	1.	Preprocess UKB data outside R if needed (e.g., ensure you have a phased VCF for APOE region, and maybe filter to relevant variants).
	2.	Run 00.EVO2_scoring_UKB_VCF_to_haplotype.R to extract haplotype info. This script might produce intermediate files listing haplotypes. Monitor its memory usage; you might need to sample or chunk if running on all 500k individuals.
	3.	Run 01.EVO2_scoring_UKB_haplotype_to_sequence.R to generate sequences and score haplotypes. You might opt to only score unique haplotypes that occur above a certain frequency to manage computational load. After this, you will have haplotype scores and each individual’s leading score.
	4.	With the scoring results in hand, you can run the analysis scripts in any order, as they are somewhat independent:
	•	If interested in AD outcome specifically, run 1.TableS10.R first to get the logistic regression result for AD (this is a single analysis yielding a table).
	•	For imaging, run 2.Fig5A_FAST.R (analysis) then 3.Fig5A_FAST_plot.R to generate the brain figure for structural MRI, and similarly 4.Fig5A_SWI.R then 5.Fig5A_SWI_plot.R for SWI. Ensure you have the MRI phenotype files (possibly one for volumes, one for SWI metrics) and that they are merged with the genetic data by participant ID in these scripts.
	•	For proteomics, run 6.Fig5B_Olink.R to perform all the protein associations. This will produce a data frame of proteins with p-values. Then run 7.Fig5BC_Olink_plot.R to create the volcano plot and perform GO/KEGG enrichment. You will need the list of significant hits; the script likely handles filtering (e.g., FDR < 0.05) internally. Ensure that the gprofiler2 queries in the script have internet access or you have a local GO database, as enrichment analysis often queries an external service.
	•	For PheWAS, run 8.Fig5D_PheWAS.R. This may take some time as it involves many regressions. You might want to run it on a computing cluster or with parallelization if possible. It will output a file or object with all results. Then run 9.Fig5D_PheWAS_plot.R to generate the Manhattan-style plot. You may need to adjust the plotting if there are extremely many phenotypes; typically one might filter to a subset or plot by categories to keep it readable.
	5.	Check all outputs: The scripts may print summaries to console and save images as PDF/PNG. Make sure to direct these outputs to a directory of your choice. After running, you should obtain the components of Figure 5:
	•	Fig 5A: brain region maps for structural (and possibly SWI) associations.
	•	Fig 5B: volcano plot of proteomic associations.
	•	Fig 5C: bar plot of pathway enrichment from proteomics.
	•	Fig 5D: Manhattan plot of PheWAS results.

These comprehensive analyses demonstrate the power of integrating the DNA foundation model’s scores with a rich dataset like UKB, revealing that the APOE locus (as characterized by the model) influences a broad swath of phenotypes, many of which align with known biology of APOE (e.g., neurodegeneration, cardiovascular traits).

Contact

For any questions, clarifications, or issues related to this code or the study, please contact:
	•	Rui Zhu – email: rui.zhu.rz399@yale.edu
	•	Xiaopu Zhou – email: xiaopu.zhou@sickkids.ca

We welcome feedback and collaboration. If you use this code or the Evo2 scoring approach in your own research, please cite our paper (once available) and acknowledge the source. Thank you!
