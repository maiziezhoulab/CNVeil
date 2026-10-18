# CNVeil benchmark results

Raw tool outputs and processed benchmark results are organized into seven independently extractable archives. Each dataset includes all available tool outputs and corresponding processed results.

## Downloads

| Archive | Datasets | Download |
|---|---|---|
| `benchmark_01.tar.gz` | S0_section_A, S0_section_B | [Zenodo](https://doi.org/10.5281/zenodo.23458030) |
| `benchmark_02.tar.gz` | S0_section_C, S0_section_D | [Zenodo](https://doi.org/10.5281/zenodo.23458042) |
| `benchmark_03.tar.gz` | S0_section_E, KTN302, T10, T16, wellDR-seq | [Zenodo](https://doi.org/10.5281/zenodo.23458049) |
| `benchmark_04.tar.gz` | TN1, TN2, p1, p3, p4, p5, shared files | [Zenodo](https://doi.org/10.5281/zenodo.23458052) |
| `benchmark_05.tar.gz` | TN3, TN4 | [Zenodo](https://doi.org/10.5281/zenodo.23458056) |
| `benchmark_06.tar.gz` | TN5, TN6 | [Zenodo](https://doi.org/10.5281/zenodo.23458064) |
| `benchmark_07.tar.gz` | TN7, TN8 | [Zenodo](https://doi.org/10.5281/zenodo.23458070) |

Each archive is limited to 50 GB (50,000,000,000 bytes).

## Tool availability

The following tool sets describe the available **raw outputs**:

- **Set A:** CNRein, CNVeil, SEACON, SPRINTER, Alleloscope, AneuFinder, CHISEL-Eagle, FLCNA, Ginkgo, HMMcopy, SCOPE, SeCNV.
- **Set B:** All tools in Set A, plus RCCAE.
- **Set C:** CNVeil, SPRINTER, AneuFinder, FLCNA, HMMcopy, RCCAE, SCOPE, SeCNV. Ginkgo simulation outputs are provided separately in the shared workbooks within archive 04.
- **Set D:** CNVeil and Ginkgo.

| Dataset | Archive | Raw tool set | Processed results |
|---|---|---|---|
| S0_section_A | 01 | A | Allelic CN and total CN |
| S0_section_B | 01 | A | Allelic CN and total CN |
| S0_section_C | 02 | A | Allelic CN and total CN |
| S0_section_D | 02 | A | Allelic CN and total CN |
| S0_section_E | 03 | A | Allelic CN and total CN |
| KTN302 | 03 | B | Total CN |
| T10 | 03 | B | Total CN |
| T16 | 03 | B | Total CN |
| wellDR-seq | 03 | D | Allelic CN and total CN |
| TN1 | 04 | A | Allelic CN and total CN |
| TN2 | 04 | A | Allelic CN and total CN |
| p1 | 04 | C | Total CN |
| p3 | 04 | C | Total CN |
| p4 | 04 | C | Total CN |
| p5 | 04 | C | Total CN |
| TN3 | 05 | A | Allelic CN and total CN |
| TN4 | 05 | A | Allelic CN and total CN |
| TN5 | 06 | A | Allelic CN and total CN |
| TN6 | 06 | A | Allelic CN and total CN |
| TN7 | 07 | A | Allelic CN and total CN |
| TN8 | 07 | A | Allelic CN and total CN |

Processed categories indicate which analyses are included; they do not imply that every raw tool has a processed result in every category.

## Result locations

Paths below are relative to the extraction directory.

| Content | Path |
|---|---|
| Raw tool outputs | `benchmark_results/<dataset>/raw/<tool>/` |
| Processed allelic-CN results | `benchmark_results/<dataset>/processed/allelicCN/result/` |
| Processed total-CN results | `benchmark_results/<dataset>/processed/total/cnv_results/` |
| Dataset overview PDFs | `benchmark_results/<dataset>/processed/total/overview/` |
| Shared files and cross-dataset summaries | `benchmark_results/shared/` |
| Per-file inventory | `benchmark_results/manifests/benchmark_NN_files.tsv` |

Tool directory names retain their original spelling:

| Tool | Directory |
|---|---|
| CNRein | `CNRein` |
| CNVeil | `CNVeil` |
| SEACON | `SEACON` |
| SPRINTER | `SPRINTER` |
| Alleloscope | `alleloscope` |
| AneuFinder | `aneufinder` |
| CHISEL-Eagle | `chisel_eagle` |
| FLCNA | `flcna` |
| Ginkgo | `ginkgo` |
| HMMcopy | `hmmcopy` |
| RCCAE | `rccae` |
| SCOPE | `scope` |
| SeCNV | `secnv` |

For example:

- TN1 raw CNVeil outputs: `benchmark_results/TN1/raw/CNVeil/`
- TN1 processed total-CN results: `benchmark_results/TN1/processed/total/cnv_results/`
- S0 section A raw Alleloscope outputs: `benchmark_results/S0_section_A/raw/alleloscope/`

Each archive contains a file inventory recording every payload file, its original source path, and its path within the archive. Consult this inventory for exact filenames.

## Ginkgo simulation workbooks

Archive 04 contains the original workbooks under:

`benchmark_results/shared/raw/ginkgo/simulation_workbooks/`

Each workbook is included once and remains unchanged.

| Dataset | Workbook | Sheet position (1-based) | Sheet name |
|---|---|---:|---|
| p1 | `sim_normal.xlsx` | 1 | Sheet1 |
| p4 | `sim_normal.xlsx` | 3 | Sheet3 |
| p3 | `sim_normal_rep.xlsx` | 1 | 3 |
| p5 | `sim_normal_rep.xlsx` | 2 | 5 |

Other workbook sheets are retained but are not assigned to a benchmark dataset here.

Cross-dataset processed summaries are also included in archive 04, under:

`benchmark_results/shared/processed/`

## Naming and selected versions

- CNVeil `sectionA`–`sectionE` correspond to `S0_section_A`–`S0_section_E`.
- RCCAE `sim1/3/4/5` correspond to `p1/3/4/5`. RCCAE TN and S0 outputs are excluded.
- SeCNV aliases: `ploidy2_1 → p1`, `mixed2_3 → p3`, `ploidy2_4 → p4`, and `mixed2_5 → p5`. Distinct source outputs are retained in separate source-named subdirectories.
- `wellDR-seq_P1` corresponds to `wellDR-seq`.
- Alleloscope section E uses the selected `S0_section_E_filteredon_visual_inspection` output.
- Processed `_same` variants and CRC2 are excluded.
- Scripts and old outputs were filtered by the original collector, except for the explicitly selected Alleloscope section E source. Existing nested archives were retained unchanged.

## Verify and extract

Place `SHA256SUMS.txt` alongside the downloaded archives, then verify their checksums:

```bash
sha256sum --ignore-missing -c SHA256SUMS.txt
```

Extract an archive:

```bash
tar -xzf benchmark_01.tar.gz
```

Each archive can be extracted independently. Extract multiple archives into the same directory to collect their datasets.
