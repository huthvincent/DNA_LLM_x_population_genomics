#!/bin/bash
# ==============================================================================
# 0.EVO2_scoring_PanGenome.sh
# Author: Xiaopu Zhou
# Description: Master pipeline to download and process human pangenome assemblies.
# Steps:
#   1. Download GFA assembly
#   2. Run minimap2 alignment to APOE reference
#   3. Concatenate PAF alignment files
#   4. Extract target region from reference
#   5. Extract aligned regions from pangenome assemblies
# ==============================================================================

# Exit on error
set -euo pipefail

# === Configuration ===
PanGenomePATH="/hpf/largeprojects/tcagstor/users/xpzhou/2025June/PanGenome"
MINIMAP2_PATH="/hpf/largeprojects/tcagstor/users/xpzhou/software/minimap2-2.26_x64-linux/minimap2"
REFERENCE="${PanGenomePATH}/minimap2/APOE_ref_44841743_44921743_GRCh38.fasta"

# === Step 1: Download GFA Assembly ===
echo ">>> Step 1: Download pangenome GFA assembly for chr19"
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr19.hprc-v1.0-pggb.gfa.gz

# === Step 2: Run minimap2 to align assemblies ===
echo ">>> Step 2: Running minimap2 alignments"

mkdir -p "${PanGenomePATH}/minimap2"

cat << 'EOF' > minimap2.sh
#!/bin/bash
set -euo pipefail

PanGenomePATH="/hpf/largeprojects/tcagstor/users/xpzhou/2025June/PanGenome"
MINIMAP2="/hpf/largeprojects/tcagstor/users/xpzhou/software/minimap2-2.26_x64-linux/minimap2"

for assembly in ${PanGenomePATH}/assmbly/*.fa.gz; do
  filename=$(basename "$assembly" .fa.gz)

  echo "Aligning: $filename"

  "$MINIMAP2" \
    -x asm20 \
    -t 20 \
    --secondary=no \
    "$assembly" \
    ${PanGenomePATH}/minimap2/APOE_ref_44841743_44921743_GRCh38.fasta \
    > ${PanGenomePATH}/minimap2/${filename}.paf
done
EOF

chmod +x minimap2.sh

# Submit job to SLURM
sbatch --mem=50G --time=144:00:00 --cpus-per-task=20 --wrap="bash minimap2.sh"

# === Step 3: Concatenate PAF files ===
echo ">>> Step 3: Concatenating minimap2 output files"

mkdir -p "${PanGenomePATH}/minimap2_concat"

cat << 'EOF' > concat_paf.sh
#!/bin/bash
set -euo pipefail

folder="/hpf/largeprojects/tcagstor/users/xpzhou/2025June/PanGenome/minimap2"
output="/hpf/largeprojects/tcagstor/users/xpzhou/2025June/PanGenome/minimap2_concat/combined_output.txt"

> "$output"

for file in "$folder"/*.paf; do
  filename=$(basename "$file")
  awk -v fname="$filename" '{print fname "\t" $0}' "$file" >> "$output"
done
EOF

chmod +x concat_paf.sh
bash concat_paf.sh

# === Step 4: Extract APOE target region from GRCh38 reference ===
echo ">>> Step 4: Extracting APOE region from GRCh38 reference"

mkdir -p "${PanGenomePATH}/APOE_locus"

samtools faidx "$REFERENCE"
samtools faidx "$REFERENCE" hg38_dna:22-79918 \
  > "${PanGenomePATH}/APOE_locus/APOE_minimap2_match_GRCh38_ref.fa"

# === Step 5: Extract APOE regions from pangenome assemblies ===
echo ">>> Step 5: Extracting APOE-mapped regions from pangenome assemblies"

cat << 'EOF' > extract.pangenome.sh
#!/bin/bash
set -euo pipefail

input_file="./APOE_locus/pangneome_extract.list"

tail -n +2 "$input_file" | while IFS=$'\t' read -r fasta_file full_region; do
  region=$(echo "$full_region")
  fasta_file=$(echo "$fasta_file" | sed -e 's/paf/fa.gz/g')

  if [[ ! -f ./assmbly/"$fasta_file" ]]; then
    echo "FASTA file not found: $fasta_file"
    continue
  fi

  echo "Extracting $region from $fasta_file"
  samtools faidx ./assmbly/"$fasta_file"
  samtools faidx ./assmbly/"$fasta_file" "$region" \
    > ./APOE_locus/"$fasta_file".fasta
done
EOF

chmod +x extract.pangenome.sh

sbatch --mem=50G --time=144:00:00 --cpus-per-task=20 --wrap="bash extract.pangenome.sh"

echo ">>> All steps submitted or completed."
