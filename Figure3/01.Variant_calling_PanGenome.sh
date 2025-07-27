#!/bin/bash
# ==============================================================================
# 01.Variant_calling_PanGenome.sh
# Author: Xiaopu Zhou
# Description: Complete APOE variant calling pipeline from alignment to VCF
# Dependencies: minimap2, samtools, bcftools
# ==============================================================================

set -euo pipefail

# === [1] Set paths ===
PanGenomePATH="/hpf/largeprojects/tcagstor/users/xpzhou/2025June/PanGenome"
MINIMAP2="/hpf/largeprojects/tcagstor/users/xpzhou/software/minimap2-2.26_x64-linux/minimap2"
REF="${PanGenomePATH}/minimap2/APOE_ref_44841743_44921743_GRCh38.fasta"
FA_DIR="${PanGenomePATH}/APOE_locus"
BAM_DIR="${PanGenomePATH}/variant_call/BAM"
BAM_SORT_DIR="${PanGenomePATH}/variant_call/BAM_sort"
VCF_DIR="${PanGenomePATH}/variant_call/VCF"

mkdir -p "$BAM_DIR" "$BAM_SORT_DIR" "$VCF_DIR"

# === [2] Align pangenome assemblies to GRCh38 APOE reference ===
echo "Step 1: Minimap2 alignment..."

for fasta in "${FA_DIR}"/*.fa.gz.fasta; do
  fname=$(basename "$fasta" .fa.gz.fasta)
  echo "Aligning ${fname}..."
  
  "$MINIMAP2" -a -x asm20 -t 20 --secondary=no "$REF" "$fasta" | \
    samtools sort -@ 10 -o "${BAM_DIR}/${fname}_GRCh38_align.bam" --write-index -

  samtools index "${BAM_DIR}/${fname}_GRCh38_align.bam"
done

# === [3] Sort BAMs (optional if needed) ===
echo "Step 2: Sorting BAMs..."

for bam in "${BAM_DIR}"/*.bam; do
  fname=$(basename "$bam" .bam)
  echo "Sorting ${fname}..."
  
  samtools sort "$bam" -o "${BAM_SORT_DIR}/${fname}.sorted.bam"
  samtools index "${BAM_SORT_DIR}/${fname}.sorted.bam"
done

# === [4] Variant calling with bcftools mpileup + call ===
echo "Step 3: Variant calling..."

# Generate list of BAMs
BAM_LIST="${VCF_DIR}/BAM.list"
find "${BAM_SORT_DIR}" -name "*.sorted.bam" > "$BAM_LIST"

# Run bcftools mpileup and call
bcftools mpileup \
  -f "$REF" \
  -b "$BAM_LIST" \
  -Ou \
  --threads 10 \
  --no-BAQ | bcftools call \
    --ploidy 1 \
    --threads 10 \
    -mv \
    -Oz \
    -o "${VCF_DIR}/1KG_Pangenome_assembly_call.APOE.vcf.gz"

# Index final VCF
bcftools index "${VCF_DIR}/1KG_Pangenome_assembly_call.APOE.vcf.gz"

echo "✅ Variant calling pipeline complete."
