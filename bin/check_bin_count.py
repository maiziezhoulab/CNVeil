import os, glob, argparse
import pandas as pd
import matplotlib.pyplot as plt
import math

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True)
parser.add_argument("-o", required=True)
args = parser.parse_args()

os.makedirs(args.o, exist_ok=True)

files = glob.glob(os.path.join(args.i, "*", "CNV_*.csv.T"))

counts = {}
for f in files:
    data = os.path.basename(os.path.dirname(f))
    s = pd.read_csv(f, usecols=[0]).iloc[:,0].astype(str).str.split(":").str[0]
    counts[data] = s.value_counts()

df = pd.DataFrame(counts).fillna(0).astype(int)

# summary
summary = pd.DataFrame({
    "min": df.min(axis=1),
    "min_tool": df.idxmin(axis=1),
    "max": df.max(axis=1),
    "max_tool": df.idxmax(axis=1),
    "median": df.median(axis=1),
})

summary["range"] = summary["max"] - summary["min"]
summary["rel_range"] = summary["range"] / summary["max"]

# save counts
df.to_csv(os.path.join(args.o, "bin_count_per_chrom.csv"))

# save summary as excel
summary.to_excel(os.path.join(args.o, "bin_count_summary.xlsx"))

# heatmap (raw counts)
plt.figure(figsize=(max(6, df.shape[1]*0.6), max(4, df.shape[0]*0.4)))
plt.imshow(df, aspect="auto")
plt.colorbar(label="Bin count")
plt.xticks(range(df.shape[1]), df.columns, rotation=90)
plt.yticks(range(df.shape[0]), df.index)
plt.title("Bin counts per chromosome")
plt.tight_layout()
plt.savefig(os.path.join(args.o, "bin_count_heatmap.pdf"))
plt.close()

# heatmap (max - current)
diff = df.max(axis=1).values[:,None] - df

plt.figure(figsize=(max(6, diff.shape[1]*0.6), max(4, diff.shape[0]*0.4)))
plt.imshow(diff, aspect="auto")
plt.colorbar(label="max - bin_count")
plt.xticks(range(diff.shape[1]), diff.columns, rotation=90)
plt.yticks(range(diff.shape[0]), diff.index)
plt.title("Missing bins relative to chromosome max")
plt.tight_layout()
plt.savefig(os.path.join(args.o, "bin_count_missing_heatmap.pdf"))
plt.close()

# histogram per chromosome
chroms = df.index
cols = 4
rows = math.ceil(len(chroms)/cols)

fig, axes = plt.subplots(rows, cols, figsize=(cols*4, rows*3))
axes = axes.flatten()

for i, chrom in enumerate(chroms):
    axes[i].hist(df.loc[chrom], bins=15)
    axes[i].set_title(chrom)

for j in range(i+1, len(axes)):
    axes[j].axis("off")

plt.tight_layout()
plt.savefig(os.path.join(args.o, "bin_count_histograms.pdf"))
plt.close()