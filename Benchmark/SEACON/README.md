# SEACON: hg19 workflow

Activate your working SEACON environment, with its SAMtools/bcftools and other
dependencies available. Then run:

```bash
conda activate SEACON
bash run.sh /path/to/hg19_cell_bams /path/to/hg19.fa \
  /path/to/hg19_mappability.bigWig /path/to/phased_normal_snps.vcf.gz \
  /path/to/hg19_matched_normal.bam /path/to/new_SEACON_TN1
```

Inputs: individual-cell BAM directory, indexed hg19 FASTA, matching mappability
bigWig, indexed phased SNP VCF, and indexed matched-normal BAM. All must use the
same actual hg19 coordinates and matching chromosome names. VCF generation and
phasing remain upstream; the commented pseudonormal workflow is not enabled.

**The original matched-normal filename contained `hg38`, while the FASTA and
mappability were hg19.** Confirm the real alignment build before supplying it.
This wrapper does not rename or convert a mismatched BAM.

| Preserved setting | Value |
|---|---|
| Read-count bins | 500,000 bp |
| Read-count processes | 12 |
| BAF block size | 50 kb |
| Matched normal | Supplied to read-count and BAF preparation |
| Component merge tolerance | 0.2 |
| Breakpoint upper filter | 5 |
| Maximum WGD | 2 |

Runs `prep_readcount`, `prep_baf`, then `call`; all other tool settings retain
installed defaults. Results and three command logs go into the new output
directory. Personal paths, inactive alternatives, and module purging are removed.
Original cluster allocation: 15 CPUs, 100 GB, 7 days.
Bash syntax checked; execution and actual assembly consistency remain unverified.

[SEACON repository](https://github.com/NabaviLab/SEACON)
