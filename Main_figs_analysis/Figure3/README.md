# eval_T10

Scripts in this folder evaluate copy-number calls on the T10 real-data benchmark.
The workflow assigns cells to T10 sectors by consensus across tools, evaluates each
tool against the T10 gold CNV profile, and plots per-sector MSE distributions.

## Files

- `t10_consensus.py`
  - Reads CNV result TSV files for `CNVeil`, `Ginkgo`, `SCOPE`, and `SeCNV`.
  - Estimates each cell's sector by comparing the cell mean copy number to fixed
    sector ploidies:
    - `A1`: 3.05
    - `A2`: 2.85
    - `H`: 1.7
    - `D`: 2.0
  - Keeps cells whose top sector assignment is supported by at least `k` tools.
    The script currently runs with `k=4`, so all four tools must agree.
  - Writes consensus outputs under `./consensus/`.
  - Takes input/output paths, tool names, ploidies, sector names, and cell-column
    name as command-line arguments.

- `t10_eval.py`
  - Reads tool CNV result files and the consensus cell-to-sector map.
  - Compares predicted per-bin copy number against the sector-specific gold CNV
    profile from `T10_gold_CNV.tsv`.
  - Computes per-cell MSE and Pearson correlation for each sector.
  - Writes per-tool metric pickles and a combined summary table under
    `./evaluation/`.
  - Takes the tool result directory, gold CNV file, consensus file, output
    directory, bin size, and cell-column name as command-line arguments.

- `t10_plot.py`
  - Reads `./evaluation/*_metrics.pkl`.
  - Computes overall weighted MSE per tool, ranks tools from lowest to highest
    MSE, and plots per-sector MSE boxplots.
  - Writes outputs under `./plot/`.
  - Takes the metrics directory, output directory, and sector order as
    command-line arguments.

## Expected Inputs

CNV result directory for consensus:

```text
/path/to/t10_cnv_results/
  CNVeil.tsv
  Ginkgo.tsv
  SCOPE.tsv
  SeCNV.tsv
```

Each tool TSV should have:

- one row per cell
- first column named `Cell`
- remaining columns named as genomic bins, for example `chr1:0-500000`
- copy-number values in the bin columns

Gold CNV file used by `t10_eval.py`, for example:

```text
/path/to/T10_gold_CNV.tsv
```

The gold file must contain:

- `CHROM`
- `start_i`
- one column per consensus sector, such as `A1`, `A2`, `H`, and `D`

## How To Run

Run these commands from the `eval_T10` folder:

```bash
cd eval_T10
```

1. Build consensus cell-to-sector assignments:

```bash
python3 t10_consensus.py \
  --input_dir /path/to/t10_cnv_results \
  --output_dir consensus \
  --tools CNVeil,Ginkgo,SCOPE,SeCNV \
  --k 4
```

Outputs:

```text
consensus/cell2sector_consensus.pkl
consensus/cell2sector_support.pkl
consensus/report.txt
```

2. Evaluate each tool against the T10 gold CNV profile:

```bash
python3 t10_eval.py \
  --input_dir /path/to/t10_cnv_results \
  --gold_tsv /path/to/T10_gold_CNV.tsv \
  --consensus_pkl consensus/cell2sector_consensus.pkl \
  --output_dir evaluation
```

Outputs:

```text
evaluation/<tool>_metrics.pkl
evaluation/summary.tsv
```

3. Plot and rank tools by weighted MSE:

```bash
python3 t10_plot.py \
  --metrics_dir evaluation \
  --output_dir plot
```

Outputs:

```text
plot/overall_weighted_mse.tsv
plot/mse_boxplots_by_sector.pdf
```

## Dependencies

The scripts use:

- Python 3
- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `matplotlib`

Install them in the active environment if needed:

```bash
pip install numpy pandas scipy scikit-learn matplotlib
```

## Notes And Caveats

- `t10_eval.py` ignores files containing `rcCAE` in the T10 result directory.
  Change this with `--exclude`, or pass `--exclude ''` to disable filtering.
- Bin starts are converted to `start_i` using a default `bin_size` of 500 kb.
  Change this with `--bin_size`.
- `t10_eval.py` expects the consensus file at
  `consensus/cell2sector_consensus.pkl`, so run `t10_consensus.py` first.
- If tool files use a different cell-column name, pass `--cell_col`.
