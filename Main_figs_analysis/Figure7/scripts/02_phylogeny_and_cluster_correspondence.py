#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

TOTAL_CLUSTER_FILE = Path("example_data/figure7/copy_number/total_cn_cluster_assignments.csv")
HAPLOTYPE_CLUSTER_FILE = Path("example_data/figure7/copy_number/haplotype_cluster_assignments.tsv")
TOTAL_TREE_FILE = Path("results/figure7/medicc2_total_cn/consensus_cn_profiles_support_tree.new")
HAPLOTYPE_TREE_FILE = Path("results/figure7/medicc2_haplotype/input_medicc2_support_tree.new")
OUTPUT_DIR = Path("results/figure7")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def read_table(path):
    with open(path) as handle:
        first = handle.readline()
    return pd.read_csv(path, sep="\t" if "\t" in first else ",")


def standardize(path):
    df = read_table(path)
    cell_col = next(c for c in ["Cell_ID", "Cell", "cell", "cell_id"] if c in df.columns)
    cluster_col = next(c for c in ["Cluster", "cluster", "cluster_id"] if c in df.columns)
    return df[[cell_col, cluster_col]].rename(columns={cell_col: "Cell_ID", cluster_col: "Cluster"})


def label_cluster(x):
    x = str(x).strip().replace("_", " ")
    if x.lower() == "diploid":
        return "diploid"
    return x if x.lower().startswith("cluster") else f"cluster {x}"


def plot_trees():
    total_tree = Phylo.read(TOTAL_TREE_FILE, "newick")
    hap_tree = Phylo.read(HAPLOTYPE_TREE_FILE, "newick")
    for tree in [total_tree, hap_tree]:
        for tip in tree.get_terminals():
            tip.name = label_cluster(tip.name)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    Phylo.draw(total_tree, axes=axes[0], do_show=False, show_confidence=False)
    axes[0].set_title("Total CN Phylogenetic Tree")
    axes[0].set_xlabel("Evolutionary Distance")
    axes[0].set_ylabel("taxa")

    Phylo.draw(hap_tree, axes=axes[1], do_show=False, show_confidence=False)
    axes[1].set_title("Haplotype-Specific CN Tree")
    axes[1].set_xlabel("Evolutionary Distance")
    axes[1].set_ylabel("taxa")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "Figure7a_phylogenetic_trees.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "Figure7a_phylogenetic_trees.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_correspondence():
    total = standardize(TOTAL_CLUSTER_FILE).rename(columns={"Cluster": "Total_cluster"})
    hap = standardize(HAPLOTYPE_CLUSTER_FILE).rename(columns={"Cluster": "Haplotype_cluster"})
    merged = total.merge(hap, on="Cell_ID")
    merged["Total_cluster"] = merged["Total_cluster"].map(label_cluster)
    merged["Haplotype_cluster"] = merged["Haplotype_cluster"].map(label_cluster)

    table = pd.crosstab(merged["Total_cluster"], merged["Haplotype_cluster"])
    table.to_csv(OUTPUT_DIR / "total_to_haplotype_cluster_correspondence.tsv", sep="\t")

    ncols = 3
    nrows = int(np.ceil(len(table) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, total_cluster in zip(axes, table.index):
        counts = table.loc[total_cluster]
        counts = counts[counts > 0]
        x = np.arange(len(counts))
        ax.bar(x, counts.values, edgecolor="black", alpha=0.75)
        ax.set_xticks(x)
        ax.set_xticklabels([v.replace("cluster ", "") for v in counts.index], rotation=45, ha="right")
        ax.set_xlabel("Haplotype Cluster")
        ax.set_ylabel("Number of Cells")
        ax.set_title(f"{total_cluster.title()} → Cell Distribution\n(n_cells={int(counts.sum())}, n_hap_clusters={len(counts)})", fontsize=10)

    for ax in axes[len(table):]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "Figure7b_cluster_correspondence.pdf", bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "Figure7b_cluster_correspondence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


plot_trees()
plot_correspondence()
