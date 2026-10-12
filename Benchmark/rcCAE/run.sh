#!/usr/bin/env bash
set -euo pipefail
# Run using the original rcCAE installation.
RCCAE_REPO="/path/to/rccae"
root_dir="/path/to/TNBC"
bam1="$root_dir/DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb"
bam2="$root_dir/DNA_bowtie2_rmdup/DNA3022_rmdup_rg_cb"
ref="/path/to/hg19.fa"
map="/path/to/hg19_mappability.bigWig"
barcodes1="$root_dir/DNA302/unique_srr_ids1.txt"
barcodes2="$root_dir/DNA302/unique_srr_ids2.txt"
output_dir="/path/to/rcCAE/KTN302"
MCR="/path/to/MCR/v91"
mkdir -p "$output_dir"
cd "$RCCAE_REPO"
log_file="$output_dir/log.txt"
: > "$log_file"
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
chrlist=$(IFS=,; echo "${chromosomes[*]}")
seconds_s=$(date +%s)
echo 'Step 1. get read counts, GC-content and mappability'
./prep/bin/prepInput -b "$bam1" -r "$ref" -m "$map" -B "$barcodes1" \
    -c "$chrlist" -s 500000 -q 0 -o "$output_dir/readcounts1.txt" >> "$log_file" 2>&1
./prep/bin/prepInput -b "$bam2" -r "$ref" -m "$map" -B "$barcodes2" \
    -c "$chrlist" -s 500000 -q 0 -o "$output_dir/readcounts2.txt" >> "$log_file" 2>&1
file1="$output_dir/readcounts1.txt"
file2="$output_dir/readcounts2.txt"
output_file="$output_dir/readcounts.txt"
head -n 1 "$file1" > "$output_file"
tail -n +2 "$file1" > "$output_file.tmp1"
tail -n +2 "$file2" | cut -d, -f4- > "$output_file.tmp2"
paste -d, "$output_file.tmp1" "$output_file.tmp2" >> "$output_file"
rm "$output_file.tmp1" "$output_file.tmp2"
echo 'Step 2. detect tumor clones and denoise read counts'
ml Anaconda3/2023.03-1
source activate rccae
python ./cae/train.py --input "$output_dir/readcounts.txt" --epochs 50 \
    --batch_size 64 --lr 0.0001 --latent_dim 3 --seed 0 \
    --output "$output_dir" >> "$log_file" 2>&1
echo 'Step 3. detect single-cell copy number alterations'
./hmm/run_SCHMM.sh "$MCR" "$output_dir/lrc.txt" "$output_dir" 10 >> "$log_file" 2>&1
seconds_e=$(date +%s)
echo "Elapsed time: $((seconds_e-seconds_s)) seconds"
