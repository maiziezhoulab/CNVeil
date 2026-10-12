#!/usr/bin/env bash
set -euo pipefail
# Ginkgo writes results in this BED directory, as in the original workflow.
dir="/path/to/hg19_cell_beds"
ginkgo="/path/to/ginkgo/cli/ginkgo.sh"
ls "$dir"/*_hg19_rmdup.bed.gz > "$dir/cell.list"
"$ginkgo" --input "$dir" --genome hg19 \
    --binning variable_500000_76_bowtie --cells "$dir/cell.list"
