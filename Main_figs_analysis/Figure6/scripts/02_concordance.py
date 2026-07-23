#!/usr/bin/env python3

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

bulk_base = "example_data/figure6/ASCAT"
sc_base = "example_data/figure6/CNVeil"
sample_info_path = "example_data/figure6/TN_sample_info.tsv"
output_root = "results/figure6"

# Figure 6:
TN_list = ["TN6"]

# For Figure S7, use:
# TN_list = ["TN1", "TN2", "TN3", "TN4", "TN5", "TN7", "TN8"]

max_cn = 8
fontsize = 14

chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

def normalize_chr_label(chrom):
    chrom = str(chrom).strip()

    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    if chrom == "23":
        chrom = "X"
    elif chrom == "24":
        chrom = "Y"

    return "chr" + chrom


def load_sample_info():
    df = pd.read_csv(sample_info_path, sep="\t")
    df.columns = df.columns.str.strip()
    return df


def get_tumor_id(sample_info, patient_id):
    row = sample_info.loc[sample_info["Patient_ID"] == patient_id]

    if row.empty:
        raise ValueError(f"{patient_id} was not found in {sample_info_path}")

    return row["Tumor_ID"].iloc[0]


def load_bulk_segments(tumor_id):
    path = Path(bulk_base) / f"{tumor_id}.segments.txt"

    df = pd.read_csv(path, sep="\t")

    df = df.rename(
        columns={
            "startpos": "start",
            "endpos": "end",
            "nMajor": "bulk_major",
            "nMinor": "bulk_minor",
        }
    )

    df["chrom"] = df["chr"].apply(normalize_chr_label)

    return df[
        ["chrom", "start", "end", "bulk_major", "bulk_minor"]
    ]


def load_pseudo_bulk(tn):
    path = Path(sc_base) / tn / "pseudo_bulk_hapCN.tsv"

    df = pd.read_csv(path, sep="\t")

    if "cn_a" in df.columns:
        df = df.rename(columns={"cn_a": "pseudo_cn_a"})

    if "cn_b" in df.columns:
        df = df.rename(columns={"cn_b": "pseudo_cn_b"})

    df["chrom"] = df["chrom"].apply(normalize_chr_label)

    return df[
        ["chrom", "start", "end", "pseudo_cn_a", "pseudo_cn_b"]
    ]

def map_bulk_to_bins(bulk_df, pseudo_df):
    rows = []

    for chrom in chrom_order:
        bins_chr = pseudo_df[pseudo_df["chrom"] == chrom]
        segs_chr = bulk_df[bulk_df["chrom"] == chrom].copy()

        if bins_chr.empty or segs_chr.empty:
            continue

        segs_chr = segs_chr.sort_values("start").reset_index(drop=True)
        segs_chr["segment_id"] = np.arange(len(segs_chr))

        for _, bin_row in bins_chr.iterrows():
            overlaps = segs_chr[
                (segs_chr["end"] > bin_row["start"])
                & (segs_chr["start"] < bin_row["end"])
            ].copy()

            if overlaps.empty:
                continue

            overlaps["overlap"] = (
                np.minimum(overlaps["end"], bin_row["end"])
                - np.maximum(overlaps["start"], bin_row["start"])
            )

            best = overlaps.sort_values("overlap", ascending=False).iloc[0]

            rows.append(
                {
                    "chrom": chrom,
                    "start": bin_row["start"],
                    "end": bin_row["end"],
                    "pseudo_cn_a": bin_row["pseudo_cn_a"],
                    "pseudo_cn_b": bin_row["pseudo_cn_b"],
                    "segment_id": best["segment_id"],
                    "bulk_major": best["bulk_major"],
                    "bulk_minor": best["bulk_minor"],
                }
            )

    return pd.DataFrame(rows)

def phase_bulk_haplotypes(mapped_df):
    df = mapped_df.copy()

    df["bulk_imbalance"] = df["bulk_major"] - df["bulk_minor"]
    df["pseudo_imbalance"] = df["pseudo_cn_a"] - df["pseudo_cn_b"]

    segment_stats = (
        df.groupby(["chrom", "segment_id"])
        .agg(
            mean_bulk_imbalance=("bulk_imbalance", "mean"),
            mean_pseudo_imbalance=("pseudo_imbalance", "mean"),
        )
        .reset_index()
    )

    segment_stats["flip"] = (
        np.sign(segment_stats["mean_bulk_imbalance"])
        != np.sign(segment_stats["mean_pseudo_imbalance"])
    )

    df = df.merge(
        segment_stats[["chrom", "segment_id", "flip"]],
        on=["chrom", "segment_id"],
        how="left",
    )

    df["bulk_allele1"] = np.where(
        df["flip"], df["bulk_minor"], df["bulk_major"]
    )

    df["bulk_allele2"] = np.where(
        df["flip"], df["bulk_major"], df["bulk_minor"]
    )

    df["pseudo_allele1"] = df["pseudo_cn_a"]
    df["pseudo_allele2"] = df["pseudo_cn_b"]

    df["bulk_total"] = df["bulk_allele1"] + df["bulk_allele2"]
    df["pseudo_total"] = df["pseudo_allele1"] + df["pseudo_allele2"]

    return df

def build_genome_coordinates(df):
    df = df.copy()

    df["chrom"] = pd.Categorical(
        df["chrom"],
        categories=chrom_order,
        ordered=True,
    )

    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

    chr_offsets = {}
    offset = 0

    for chrom in chrom_order:
        sub = df[df["chrom"] == chrom]

        if sub.empty:
            continue

        chr_offsets[chrom] = offset
        offset += sub["end"].max()

    df["x"] = df.apply(
        lambda row: chr_offsets[str(row["chrom"])] + row["start"],
        axis=1,
    )

    return df, chr_offsets


def set_chr_axis(ax, df, chr_offsets):
    positions = []
    labels = []

    for chrom in chrom_order:
        if chrom not in chr_offsets:
            continue

        sub = df[df["chrom"] == chrom]

        if sub.empty:
            continue

        start = chr_offsets[chrom]
        width = sub["end"].max()

        ax.axvline(
            start,
            color="gray",
            linestyle="--",
            alpha=0.3,
            linewidth=0.7,
        )

        positions.append(start + width / 2)
        labels.append(chrom.replace("chr", ""))

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=fontsize - 2)
    ax.tick_params(axis="x", length=0)

def safe_corr(x, y):
    valid = pd.concat([x, y], axis=1).dropna()

    if len(valid) < 2:
        return np.nan

    if valid.iloc[:, 0].nunique() <= 1:
        return np.nan

    if valid.iloc[:, 1].nunique() <= 1:
        return np.nan

    return pearsonr(valid.iloc[:, 0], valid.iloc[:, 1])[0]


def compute_per_chr_correlations(df):
    rows = []

    for chrom in chrom_order:
        sub = df[df["chrom"] == chrom]

        if len(sub) < 2:
            continue

        rows.append(
            {
                "chrom": chrom,
                "corr_total": safe_corr(
                    sub["bulk_total"],
                    sub["pseudo_total"],
                ),
                "corr_allele1": safe_corr(
                    sub["bulk_allele1"],
                    sub["pseudo_allele1"],
                ),
                "corr_allele2": safe_corr(
                    sub["bulk_allele2"],
                    sub["pseudo_allele2"],
                ),
            }
        )

    return pd.DataFrame(rows)

def plot_05_06_07(df, corr_df, tn):
    df, chr_offsets = build_genome_coordinates(df)

    fig = plt.figure(figsize=(22, 13))

    grid = fig.add_gridspec(
        4,
        1,
        height_ratios=[2, 2, 2, 2.4],
        hspace=0.45,
    )

    ax_total = fig.add_subplot(grid[0])
    ax_allele1 = fig.add_subplot(grid[1], sharex=ax_total)
    ax_allele2 = fig.add_subplot(grid[2], sharex=ax_total)
    ax_corr = fig.add_subplot(grid[3])

    # 05: total CN
    ax_total.step(
        df["x"],
        df["bulk_total"],
        where="mid",
        label="ASCAT bulk total",
        linewidth=2,
    )

    ax_total.step(
        df["x"],
        df["pseudo_total"],
        where="mid",
        label="CNVeil pseudo total",
        linewidth=2,
    )

    ax_total.set_ylabel("Total CN", fontsize=fontsize)
    ax_total.set_title(
        f"{tn}: ASCAT bulk vs CNVeil pseudo total CN",
        fontsize=fontsize + 2,
    )

    # 06: allele1
    ax_allele1.step(
        df["x"],
        df["bulk_allele1"],
        where="mid",
        label="ASCAT bulk allele1",
        linewidth=2,
    )

    ax_allele1.step(
        df["x"],
        df["pseudo_allele1"],
        where="mid",
        label="CNVeil pseudo allele1",
        linewidth=2,
    )

    ax_allele1.set_ylabel("Allele1 CN", fontsize=fontsize)
    ax_allele1.set_title(
        f"{tn}: Allele1 CN track",
        fontsize=fontsize + 2,
    )

    # 06: allele2
    ax_allele2.step(
        df["x"],
        df["bulk_allele2"],
        where="mid",
        label="ASCAT bulk allele2",
        linewidth=2,
    )

    ax_allele2.step(
        df["x"],
        df["pseudo_allele2"],
        where="mid",
        label="CNVeil pseudo allele2",
        linewidth=2,
    )

    ax_allele2.set_ylabel("Allele2 CN", fontsize=fontsize)
    ax_allele2.set_title(
        f"{tn}: Allele2 CN track",
        fontsize=fontsize + 2,
    )

    for ax in [ax_total, ax_allele1, ax_allele2]:
        ax.set_ylim(-0.5, max_cn)
        ax.set_yticks(range(0, max_cn + 1))
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            fontsize=fontsize,
        )

    ax_total.tick_params(labelbottom=False)
    ax_allele1.tick_params(labelbottom=False)

    set_chr_axis(ax_allele2, df, chr_offsets)
    ax_allele2.set_xlabel("Chromosome", fontsize=fontsize)

    # 07: chromosome correlations
    x = np.arange(len(corr_df))
    width = 0.26

    ax_corr.bar(
        x - width,
        corr_df["corr_total"],
        width,
        label="corr_total",
    )

    ax_corr.bar(
        x,
        corr_df["corr_allele1"],
        width,
        label="corr_allele1",
    )

    ax_corr.bar(
        x + width,
        corr_df["corr_allele2"],
        width,
        label="corr_allele2",
    )

    ax_corr.axhline(0, color="black", linewidth=1)
    ax_corr.set_ylim(-0.1, 1.05)

    ax_corr.set_xticks(x)
    ax_corr.set_xticklabels(
        corr_df["chrom"].str.replace("chr", "", regex=False),
        rotation=90,
    )

    ax_corr.set_xlabel("Chromosome", fontsize=fontsize)
    ax_corr.set_ylabel("Pearson r", fontsize=fontsize)
    ax_corr.set_title(
        f"{tn}: per-chromosome correlations",
        fontsize=fontsize + 2,
    )

    ax_corr.legend(
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        fontsize=fontsize,
    )

    output_dir = Path(output_root) / tn
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.subplots_adjust(right=0.82)

    fig.savefig(
        output_dir / f"{tn}_05_06_07_stacked.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        output_dir / f"{tn}_05_06_07_stacked.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    corr_df.to_csv(
        output_dir / f"{tn}_per_chromosome_correlations.tsv",
        sep="\t",
        index=False,
    )

    print(f"Finished {tn}")

sample_info = load_sample_info()

for tn in TN_list:
    print(f"Processing {tn}")

    tumor_id = get_tumor_id(sample_info, tn)

    bulk_df = load_bulk_segments(tumor_id)
    pseudo_df = load_pseudo_bulk(tn)

    mapped_df = map_bulk_to_bins(bulk_df, pseudo_df)

    if mapped_df.empty:
        print(f"No overlapping bins for {tn}")
        continue

    phased_df = phase_bulk_haplotypes(mapped_df)
    corr_df = compute_per_chr_correlations(phased_df)

    plot_05_06_07(
        phased_df,
        corr_df,
        tn,
    )
