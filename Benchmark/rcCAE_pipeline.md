root_dir="/data/maiziezhou_lab/Datasets/singlecell_data/TNBC/"
bam1=$root_dir"DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb"
bam2=$root_dir"DNA_bowtie2_rmdup/DNA3022_rmdup_rg_cb"
bam=($root_dir"DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb" $root_dir"DNA_bowtie2_rmdup/DNA3022_rmdup_rg_cb")
ref="/data/maiziezhou_lab/yikang/SeCNV/Scripts/hg38.fa"
map="/data/maiziezhou_lab/yikang/SeCNV/Scripts/hg38_mappability.bigWig"
barcodes1=$root_dir"/DNA302/unique_srr_ids1.txt"
barcodes2=$root_dir"/DNA302/unique_srr_ids2.txt"
barcodes=$root_dir"/DNA302/unique_srr_ids.txt"
output_dir="/data/maiziezhou_lab/Weiman/single_cell_project/rccae/KTN302"log_file=$output_dir/log.txt

chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

# Convert the array to a comma-separated string
chrlist=$(IFS=,; echo "${chromosomes[*]}")

# chroms=1,2,3

test -e $log_file && rm $log_file

test ! -e $output_dir && mkdir -p $output_dir

current=`date "+%Y-%m-%d %H:%M:%S"`
seconds_s=`date -d "$current" +%s`
echo "Step 1. get read counts, GC-content and mappability......"
./prep/bin/prepInput -b $bam1 -r $ref -m $map -B $barcodes1 -c $chrlist -s 500000 -q 0 -o $output_dir/readcounts1.txt > $log_file 2>&1
./prep/bin/prepInput -b $bam2 -r $ref -m $map -B $barcodes2 -c $chrlist -s 500000 -q 0 -o $output_dir/readcounts2.txt > $log_file 2>&1
file1=$output_dir/readcounts1.txt
file2=$output_dir/readcounts2.txt
output_file=$output_dir/readcounts.txt
header=$(head -n 1 "$file1")
echo "$header" > "$output_file"
tail -n +2 "$file1" > "$output_file.tmp1"
tail -n +2 "$file2" | cut -d, -f4- > "$output_file.tmp2"
paste -d, "$output_file.tmp1" "$output_file.tmp2" >> "$output_file"
rm "$output_file.tmp1" "$output_file.tmp2"

echo "Step 2. detect tumor clones and denoise read counts data......"
ml Anaconda3/2023.03-1
# module load CUDA/11.7.0
source activate rccae
python ./cae/train.py --input $output_dir/readcounts.txt --epochs 50 --batch_size 64 --lr 0.0001 --latent_dim 3 --seed 0 --output $output_dir >> $log_file 2>&1

echo "Step 3. detect single-cell copy number alterations......"
# the path to MCR v91 needs to be specified
# module load MATLAB/2023b
./hmm/run_SCHMM.sh /data/maiziezhou_lab/Weiman/software/MCR/v91 $output_dir/lrc.txt $output_dir 10 >> $log_file 2>&1

current=`date "+%Y-%m-%d %H:%M:%S"`
seconds_e=`date -d "$current" +%s`
let secs=seconds_e-seconds_s
echo "Elapsed time: $secs seconds!"

exit 0

