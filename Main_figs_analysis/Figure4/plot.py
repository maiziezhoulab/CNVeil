import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import seaborn as sns
from math import pi
import os
import warnings
warnings.filterwarnings('ignore')
  
INPUT_FILE = '/data/maiziezhou_lab/Weiman/CNVeil/analysis/CNVeil_scDNA-seq process_06_06_2025.xlsx'
 
OUTPUT_DIR = './figures'
 
os.makedirs(OUTPUT_DIR, exist_ok=True)


plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['figure.dpi'] = 150
 
TOOL_COLORS = {
    'FACS': '#2d3436',
    '10X infer': '#2d3436',
    'CNVeil': '#e74c3c',
    # 'FLCNA': '#c0392b',
    'CNVeil_AF_FEB': '#e74c3c',
    'CNVeil_AF_July': '#c0392b',
    'AF_reduce': '#e67e22',
    'Ginkgo': '#3498db',
    'AneuFinder': '#9b59b6',
    'FLCNA': '#1abc9c',
    'CHISEL_whatshap': '#f39c12',
    'CHISEL_eagle': '#d35400',
    'CHISEL': '#d35400',
    'HMMcopy': '#27ae60',
    'SPRINTER': '#e91e63',
    'Alleloscope': '#00bcd4',
    'SeCNV': '#795548',
    'CNRein': '#607d8b',
    'SEACON': '#ff5722',
    'SCOPE': '#673ab7',
}

def get_tool_color(tool):
    """Get color for a tool, with fallback"""
    return TOOL_COLORS.get(tool, '#95a5a6')
 

xlsx_path = INPUT_FILE 
 
# Read Excel TN sheet ONLY to get FACS (ground truth)
df_tn_excel = pd.read_excel(xlsx_path, sheet_name='Big Data Ploidy')
df_tn_excel = df_tn_excel.set_index('Dataset')

# Keep only FACS if present
tn_ref_cols = []
if 'FACS' in df_tn_excel.columns:
    tn_ref_cols.append('FACS')
df_tn_ref = df_tn_excel[tn_ref_cols].copy() if tn_ref_cols else pd.DataFrame(index=df_tn_excel.index)

# Read TSV summary and pivot to wide matrix: Dataset x Tool
df_ploidy = pd.read_csv(INPUT_FILE, sep='\t')

# Standardize column names just in case
expected_cols = {'Dataset', 'Tool', 'MeanPloidy'}
missing = expected_cols - set(df_ploidy.columns)
if missing:
    raise ValueError(f"Ploidy TSV missing columns: {missing}. Found: {df_ploidy.columns.tolist()}")

df_ploidy['MeanPloidy'] = pd.to_numeric(df_ploidy['MeanPloidy'], errors='coerce')

# Pivot: Dataset as index, tools as columns
df_tn_tools = df_ploidy.pivot_table(
    index='Dataset',
    columns='Tool',
    values='MeanPloidy',
    aggfunc='mean'
)

# Merge FACS (Excel) + tools (TSV)
df_tn = df_tn_ref.join(df_tn_tools, how='outer')

# Convert everything to numeric
for col in df_tn.columns:
    df_tn[col] = pd.to_numeric(df_tn[col], errors='coerce')
 
df_10x = pd.read_excel(xlsx_path, sheet_name='Sheet1')
df_10x = df_10x.rename(columns={'10X Dataset': 'Dataset'})
df_10x = df_10x.loc[:, ~df_10x.columns.str.contains('Unnamed')]
df_10x = df_10x.set_index('Dataset')
for col in df_10x.columns:
    df_10x[col] = pd.to_numeric(df_10x[col], errors='coerce')

# -----------------------------
# aCGH metrics (keep original way from Excel)
# -----------------------------
df_acgh_raw = pd.read_excel(xlsx_path, sheet_name='T10 acgh')
 
def parse_acgh_subclone(df):
    records = []
    current_metric = None
    current_subclone = None

    for idx, row in df.iterrows():
        if pd.notna(row['by subclone']) and row['by subclone'] != 'Metric':
            current_metric = row['by subclone']
        if pd.notna(row['Unnamed: 1']) and row['Unnamed: 1'] != 'Subclone':
            current_subclone = row['Unnamed: 1']
        if pd.notna(row['Unnamed: 2']) and pd.notna(row['Unnamed: 3']):
            try:
                value = float(row['Unnamed: 2'])
                tool = row['Unnamed: 3']
                if current_metric and current_subclone:
                    records.append({
                        'Metric': current_metric,
                        'Subclone': current_subclone,
                        'Value': value,
                        'Tool': tool
                    })
            except:
                pass
    return pd.DataFrame(records)

# Parse the right section (overall)
def parse_acgh_overall(df):
    records = []
    current_metric = None

    for idx, row in df.iterrows():
        if pd.notna(row['overall']) and row['overall'] != 'Metric':
            current_metric = row['overall']
        if pd.notna(row['Unnamed: 6']) and pd.notna(row['Unnamed: 7']):
            try:
                value = float(row['Unnamed: 6'])
                tool = row['Unnamed: 7']
                if current_metric:
                    records.append({
                        'Metric': current_metric,
                        'Value': value,
                        'Tool': tool
                    })
            except:
                pass
    return pd.DataFrame(records)

df_subclone = parse_acgh_subclone(df_acgh_raw)
df_overall = parse_acgh_overall(df_acgh_raw)

print(f"TN samples shape: {df_tn.shape}")
print(f"10X samples shape: {df_10x.shape}")
print(f"Subclone metrics: {len(df_subclone)} records")
print(f"Overall metrics: {len(df_overall)} records") 
 
def plot_heatmap():
    """Create heatmap of ploidy estimates (TN1-8 and 10X only) with center=2 (white)."""
    import matplotlib.colors as mcolors
    from matplotlib.colors import LinearSegmentedColormap

    tn_wanted = [f"TN{i}" for i in range(1, 9)]
    tn_idx = [x for x in tn_wanted if x in df_tn.index]
    if len(tn_idx) == 0:
        raise ValueError("None of TN1..TN8 found in df_tn.index. Check dataset names in Big Data Ploidy / TSV.")

    # 10X: keep all rows in df_10x (assumes df_10x contains only 10X datasets already)
    if df_10x.shape[0] == 0:
        raise ValueError("df_10x is empty. Check your Sheet1 reading logic.")

    
    tn_tools = ['FACS', 'CNVeil', 'AneuFinder', 'SPRINTER', 'Ginkgo', 'CHISEL', 
                'SeCNV', 'Alleloscope', 'FLCNA','HMMcopy', 'CNRein', 'SEACON']
    x10_tools = ['10X infer', 'Ginkgo', 'AneuFinder', 'FLCNA', 'HMMcopy',
                 'CNVeil_AF_July', 'SPRINTER', 'SeCNV', 'SEACON', 'AF_reduce', 'SCOPE']

    tn_subset = df_tn.loc[tn_idx, [c for c in tn_tools if c in df_tn.columns]]
    x10_subset = df_10x.loc[:, [c for c in x10_tools if c in df_10x.columns]]

    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        ["#2166ac", "#ffffff", "#b2182b"],  # blue, white, red (classic diverging)
        N=256
    )
 
    tn_arr  = tn_subset.to_numpy(dtype=float)
    x10_arr = x10_subset.to_numpy(dtype=float)

    tn_min  = np.nanmin(tn_arr)  if np.isfinite(np.nanmin(tn_arr))  else np.nan
    tn_max  = np.nanmax(tn_arr)  if np.isfinite(np.nanmax(tn_arr))  else np.nan
    x10_min = np.nanmin(x10_arr) if np.isfinite(np.nanmin(x10_arr)) else np.nan
    x10_max = np.nanmax(x10_arr) if np.isfinite(np.nanmax(x10_arr)) else np.nan

    global_min = np.nanmin([tn_min, x10_min])
    global_max = np.nanmax([tn_max, x10_max])

    # Ensure vmin <= 0-ish and vcenter=2 always works
    vmin = min(0.0, float(global_min)) if np.isfinite(global_min) else 0.0
    vmax = float(global_max) if np.isfinite(global_max) else 6.0
    if vmax <= 2.0:
        vmax = 2.1  # avoid degenerate norm

    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=2.0, vmax=vmax)
 
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # TN heatmap
    ax1 = axes[0]
    sns.heatmap(
        tn_subset,
        annot=True, fmt='.2f',
        cmap=cmap, norm=norm,
        ax=ax1,
        cbar_kws={'label': 'Ploidy (center=2)', 'shrink': 0.8},
        linewidths=0.5, linecolor='white'
    )
    ax1.set_title('TN1–TN8: Ploidy Estimates by Tool', fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel('Tool', fontsize=11)
    ax1.set_ylabel('Dataset', fontsize=11)
    ax1.tick_params(axis='x', rotation=45)

    # 10X heatmap
    ax2 = axes[1]
    sns.heatmap(
        x10_subset,
        annot=True, fmt='.2f',
        cmap=cmap, norm=norm,
        ax=ax2,
        cbar_kws={'label': 'Ploidy (center=2)', 'shrink': 0.8},
        linewidths=0.5, linecolor='white'
    )
    ax2.set_title('10X: Ploidy Estimates by Tool', fontsize=14, fontweight='bold', pad=15)
    ax2.set_xlabel('Tool', fontsize=11)
    ax2.set_ylabel('Dataset', fontsize=11)
    ax2.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig1_heatmap.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f'{OUTPUT_DIR}/fig1_heatmap.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Figure 1: Heatmap (TN1-8 + 10X) saved with center=2 colormap")
  
def plot_barplot():
    """Create bar plot comparing tool errors from ground truth"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # TN samples - error from FACS
    ax1 = axes[0]
    tools_tn = [c for c in df_tn.columns if c != 'FACS']
    errors_tn = {}
    for tool in tools_tn:
        valid = df_tn[[tool, 'FACS']].dropna()
        if len(valid) > 0:
            mae = np.mean(np.abs(valid[tool] - valid['FACS']))
            errors_tn[tool] = mae
    
    errors_tn = dict(sorted(errors_tn.items(), key=lambda x: x[1]))
    colors = [get_tool_color(t) for t in errors_tn.keys()]
    bars1 = ax1.barh(list(errors_tn.keys()), list(errors_tn.values()), color=colors, edgecolor='white')
    ax1.set_xlabel('Mean Absolute Error (vs FACS)', fontsize=11)
    ax1.set_title('TN Samples: Tool Accuracy\n(lower is better)', fontsize=14, fontweight='bold')
    ax1.axvline(x=0, color='black', linewidth=0.5)
    
    # Add value labels
    for bar, val in zip(bars1, errors_tn.values()):
        ax1.text(val + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.3f}', 
                va='center', fontsize=9)
    
    # 10X samples - error from 10X infer
    ax2 = axes[1]
    tools_10x = [c for c in df_10x.columns if c != '10X infer']
    errors_10x = {}
    for tool in tools_10x:
        valid = df_10x[[tool, '10X infer']].dropna()
        if len(valid) > 0:
            mae = np.mean(np.abs(valid[tool] - valid['10X infer']))
            errors_10x[tool] = mae
    
    errors_10x = dict(sorted(errors_10x.items(), key=lambda x: x[1]))
    colors = [get_tool_color(t) for t in errors_10x.keys()]
    bars2 = ax2.barh(list(errors_10x.keys()), list(errors_10x.values()), color=colors, edgecolor='white')
    ax2.set_xlabel('Mean Absolute Error (vs 10X infer)', fontsize=11)
    ax2.set_title('10X Samples: Tool Accuracy\n(lower is better)', fontsize=14, fontweight='bold')
    ax2.axvline(x=0, color='black', linewidth=0.5)
    
    for bar, val in zip(bars2, errors_10x.values()):
        ax2.text(val + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.3f}',
                va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig2_barplot_error.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f'{OUTPUT_DIR}/fig2_barplot_error.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Figure 2: Bar plot saved")

def plot_scatter():
    """Scatter plots of predicted vs FACS ploidy, panels sorted by MAE (low to high)."""
    key_tools_tn = [
        'CNVeil', 'AneuFinder', 'SPRINTER', 'Ginkgo', 'CHISEL',
        'SeCNV', 'Alleloscope', 'FLCNA', 'HMMcopy', 'CNRein'
    ]

    stats = {}
    for tool in key_tools_tn:
        if tool in df_tn.columns:
            valid = df_tn[['FACS', tool]].dropna()
            if len(valid) > 0:
                r = valid['FACS'].corr(valid[tool])
                mae = np.mean(np.abs(valid['FACS'] - valid[tool]))
                stats[tool] = {'r': r, 'mae': mae}

    sorted_tools = sorted(stats.keys(), key=lambda t: stats[t]['mae'])

    fig, axes = plt.subplots(2, 5, figsize=(16, 10))

    for idx, tool in enumerate(sorted_tools):
        ax = axes[idx // 5, idx % 5]
        valid = df_tn[['FACS', tool]].dropna()
        ax.scatter(valid['FACS'], valid[tool], c=get_tool_color(tool),
                   s=80, alpha=0.8, edgecolor='white', linewidth=1)

        lims = [min(valid['FACS'].min(), valid[tool].min()) - 0.2,
                max(valid['FACS'].max(), valid[tool].max()) + 0.2]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1, label='y=x')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        r = stats[tool]['r']
        mae = stats[tool]['mae']
        # MAE first, then r
        ax.set_title(f'{tool}\nMAE={mae:.3f}, r={r:.3f}', fontsize=11)
        ax.set_xlabel('FACS Ploidy', fontsize=10)
        ax.set_ylabel(f'{tool} Ploidy', fontsize=10)

        for sample in valid.index:
            ax.annotate(sample, (valid.loc[sample, 'FACS'], valid.loc[sample, tool]),
                        fontsize=7, alpha=0.7, xytext=(3, 3), textcoords='offset points')

    for idx in range(len(sorted_tools), 10):
        axes[idx // 5, idx % 5].axis('off')

    plt.suptitle('TN Samples: Tool Predictions vs FACS Ground Truth',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig3_scatter.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f'{OUTPUT_DIR}/fig3_scatter.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Figure 3: Scatter plot saved (sorted by MAE)")

if __name__ == "__main__":
    print("\n" + "="*60)
    print("Generating CNVeil Benchmark Figures")
    print("="*60 + "\n")
    
    plot_heatmap()
    plot_barplot()
    plot_scatter() 
    
    print("\n" + "="*60)
    print("All figures generated successfully!")
    print("="*60)
    print(f"\nOutput directory: {os.path.abspath(OUTPUT_DIR)}")
    print("\nGenerated files:")
    print("  - fig1_heatmap.png/pdf")
    print("  - fig2_barplot_error.png/pdf")
    print("  - fig3_scatter.png/pdf") 
