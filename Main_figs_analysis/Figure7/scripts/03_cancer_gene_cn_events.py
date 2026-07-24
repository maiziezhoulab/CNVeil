Library
/
03_cancer_gene_cn_events.py


#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CN_FILE = Path("example_data/figure7/copy_number/input_medicc2.tsv")
GENE_ANNOTATION_FILE = Path("example_data/figure7/gencode_v19_gene_pos.symbol.txt")
OUTPUT_DIR = Path("results/figure7")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CANCER_GENES = [
    "AURKA", "BRCA1", "CDH1", "CHEK2", "ESR1", "FHIT", "IGF1R",
    "MAP3K1", "MDM2", "MTDH", "NF2", "PGR", "PIK3CA", "RB1", "TP53"
]

cn = pd.read_csv(CN_FILE, sep="\t")
genes = pd.read_csv(GENE_ANNOTATION_FILE, sep="\t", header=None, names=["gene", "chrom", "start", "end"])


def cluster_id(x):
    x = str(x).strip().lower().replace(" ", "_")
    if x == "diploid":
        return np.nan
    for prefix in ["cluster_", "sample_"]:
        if x.startswith(prefix):
            x = x[len(prefix):]
    try:
        return int(x)
    except ValueError:
        return np.nan


cn["cluster_id"] = cn["sample_id"].map(cluster_id)
cn = cn[cn["cluster_id"].notna()].copy()
cn["cluster_id"] = cn["cluster_id"].astype(int)
cn["total_cn"] = cn["cn_a"] + cn["cn_b"]
clusters = sorted(cn["cluster_id"].unique())

records = []
for chrom, gene_set in genes[genes["gene"].isin(CANCER_GENES)].groupby("chrom"):
    chrom_cn = cn[cn["chrom"] == chrom]
    for _, gene in gene_set.iterrows():
        overlap = chrom_cn[(chrom_cn["start"] < gene["end"]) & (chrom_cn["end"] > gene["start"])]
        for cid, sub in overlap.groupby("cluster_id"):
            records.append({
                "gene": gene["gene"],
                "cluster_id": cid,
                "cn_a": sub["cn_a"].median(),
                "cn_b": sub["cn_b"].median(),
                "total_cn": sub["total_cn"].median(),
            })

gene_cn = pd.DataFrame(records)

heatmap = gene_cn.pivot_table(index="gene", columns="cluster_id", values="total_cn", aggfunc="median")
heatmap = heatmap.reindex(index=CANCER_GENES, columns=clusters).fillna(2)

fig, ax = plt.subplots(figsize=(max(8, len(clusters) * 0.55), max(5, len(CANCER_GENES) * 0.38)))
image = ax.imshow(heatmap.values, aspect="auto", cmap="OrRd", vmin=0, vmax=max(6, np.nanmax(heatmap.values)))
ax.set_xticks(np.arange(len(heatmap.columns)))
ax.set_xticklabels(heatmap.columns)
ax.set_yticks(np.arange(len(heatmap.index)))
ax.set_yticklabels(heatmap.index)
ax.set_xlabel("Cluster ID")
ax.set_ylabel("Cancer Gene")
ax.set_title("Cancer Gene Copy Number per Cluster\n(diploid = 2; higher values indicate gain)")
for i in range(heatmap.shape[0]):
    for j in range(heatmap.shape[1]):
        ax.text(j, i, f"{heatmap.iloc[i, j]:.0f}", ha="center", va="center", fontsize=7)
fig.colorbar(image, ax=ax, label="Median Total CN (cn_a + cn_b)")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "Figure7c_cancer_gene_cn_heatmap.pdf", bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "Figure7c_cancer_gene_cn_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig)

gene_cn["Amplification"] = gene_cn["total_cn"] > 2
gene_cn["LOH"] = (gene_cn["cn_a"] == 0) | (gene_cn["cn_b"] == 0)
gene_cn["CN-LOH"] = gene_cn["LOH"] & (gene_cn["total_cn"] == 2)
frequency = gene_cn.groupby("gene")[["Amplification", "CN-LOH", "LOH"]].mean().reindex(CANCER_GENES).fillna(0)
frequency.to_csv(OUTPUT_DIR / "cancer_gene_event_frequency.tsv", sep="\t")

plot_genes = frequency.index[::-1]
y = np.arange(len(plot_genes))
h = 0.24
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(y + h, frequency.loc[plot_genes, "Amplification"], h, label="Amplification")
ax.barh(y, frequency.loc[plot_genes, "CN-LOH"], h, label="CN-LOH")
ax.barh(y - h, frequency.loc[plot_genes, "LOH"], h, label="LOH")
ax.set_yticks(y)
ax.set_yticklabels(plot_genes)
ax.set_xlim(0, 1)
ax.set_xlabel("Fraction of Clusters")
ax.set_title("CNV Event Types per Cancer Gene")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "Figure7d_cancer_gene_event_frequency.pdf", bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "Figure7d_cancer_gene_event_frequency.png", dpi=300, bbox_inches="tight")
plt.close(fig)
