import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import pickle
import glob
import os 
# --- Sort cells (rows) by their average CNV ---
def plot(cnv, file=None, sort= False):
    cnv = np.clip(cnv, None, 7)
    new_row = np.full((1, cnv.shape[1]), -3)
    
    cnv = np.vstack([cnv, new_row])

    if sort:
        row_means = cnv.mean(axis=1)
        sorted_idx = np.argsort(row_means)
        cnv_sorted = cnv[sorted_idx, :]
    else:
        cnv_sorted = cnv
    print("start converting")
    cnv_sorted = cnv_sorted.astype(float)
    print("\n===========\n")
    print("type:", type(cnv_sorted))
    print("dtype:", cnv_sorted.dtype)
    print("shape:", cnv_sorted.shape)

    # Plot heatmap
    plt.figure(figsize=(12,6))
    plt.imshow(cnv_sorted, aspect='auto', cmap='RdBu_r')
    plt.colorbar(label="CNV value")
    plt.xlabel("Bins")
    if sort:
        plt.ylabel("Cells (sorted by avg CNV)")
        plt.title("CNV Matrix Heatmap (Cells Sorted by Average CNV)")
    else:
        plt.ylabel("Cells (original order)")
        plt.title("CNV Matrix Heatmap")
    if file is not None:
        plt.savefig(file)
    # plt.show()
    plt.close()
    



def smooth_total_cn(cnveil_out_dir, ref_type):
    cnv_file = glob.glob(f"{cnveil_out_dir}/CNV_*.csv")
    sample = cnv_file.split('_')[1].split('.')[0]

    print('sample: ',sample)
    # if sample.startswith('s'):
    #     id = sample[-1]
    #     sample_in = "S0_section_"+id
    # else:
    #     sample_in = sample
    # df= pd.read_csv(f"/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/Pipeline_2026_3_13/{sample_in}/CNV_{sample_in}.csv")
    #------------smooth by 5M
    df = pd.read_csv(cnv_file)
    # # Transpose
    df_t = df.T
    cells = [x.split('-')[0] for  x in df_t.iloc[0] ]
    # Make the first row the new column names
    df_t.columns = cells # set first row as header
    df_t = df_t.drop(df_t.index[0])  # drop that row from the data
    df_t['chrom'] = df_t.index.map(lambda x: x.split(':')[0])
    df_t['i'] = df_t.index.map(lambda x: (int(x.split('-')[1]) - 1 )//5000000)

    # Load chromosome sizes
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    chrom_sizes = pd.read_csv(f"{SCRIPT_DIR}/hg{ref_type}_ref_files/hg{ref_type}.chrom.sizes", sep="\t", header=None, names=["chrom", "length"])

    # Dictionary
    chr_lengths = dict(zip(chrom_sizes["chrom"], chrom_sizes["length"]))


    # Bin size
    bin_size = 5_000_000

    # Compute start and end
    df_t["start"] = df_t["i"] * bin_size
    df_t["end"] = df_t.apply(
        lambda row: min((row["i"] + 1) * bin_size, chr_lengths[row["chrom"]]),
        axis=1
    )
    # convert CNV columns to numeric
    # df_t[cells] = df_t[cells].apply(pd.to_numeric, errors="coerce")
    # df_t[cells] = df_t[cells].astype(float)
    # print(df_t.head())
    # print(df_t.iloc[0,0],type(df_t.iloc[0,0]))
    df_t_5m = df_t.groupby(['chrom','start','end'])[cells].mean()

    # Round to nearest integer
    df_t_5m = df_t_5m.round()

    # Extract numeric part for sorting, put X/Y/MT at the end
    def chrom_key(c):
        c = c.replace("chr", "")
        if c.isdigit():
            return (0, int(c))      # numeric chroms
        elif c == "X":
            return (1, 23)
        elif c == "Y":
            return (1, 24)
        elif c in ["M", "MT"]:
            return (1, 25)
        else:
            return (2, c)           # anything else

    # Sort
    df_sorted = df_t_5m.sort_values(
        by=["chrom", "start", "end"],
        key=lambda col: col.map(chrom_key) if col.name == "chrom" else col
    )
    df_sorted.to_csv(f"{cnveil_out_dir}/{sample}_smooth_5M.csv", index = False)
    plot(df_sorted[cells].values.T, f'{cnveil_out_dir}/{sample}_smooth_5M.pdf')

    #------------smooth by seg
    cnv = df_sorted.values.T
    bnd_sup_list = (1-((cnv[:,:-1] - cnv[:,1:])==0)).sum(0)
    # plt.hist(bnd_sup_list, bins = 50)
    # plt.show()
    n_tumor = (cnv.mean(1)>2.2).sum()
    n = int(n_tumor * 0.1)
    print(f"# tumor cells: {n_tumor}; # bnd support threshold: {n}")
    print("num of raw segments:",(bnd_sup_list>0).sum()+1)
    print("num of segments after filter:",(bnd_sup_list>n).sum()+1)
    bnds = [0] + list(np.where(bnd_sup_list>n)[0]+1) + [cnv.shape[1]]
    print("convert again")
    cnv = cnv.astype(float)
    print("convert finished")
    # Copy original matrix
    cnv_new = cnv.copy()

    # Replace bins within each BND with the interval's average per cell
    for i in range(len(bnds)-1):
        start = bnds[i]
        end = bnds[i+1]
        avg_values = cnv[:, start:end].mean(axis=1).round()
        cnv_new[:, start:end] = avg_values[:, None].repeat(end-start, axis=1)


    # plot(cnv_new,f'section{sample}_smooth_seg.pdf')
    plot(cnv_new,f'{cnveil_out_dir}/{sample}_smooth_seg.pdf')

    df1 = df_sorted.copy()
    df1.iloc[:,:] = cnv_new.T
    df1= df1.reset_index()

    df1.to_csv(f"{cnveil_out_dir}/{sample}_smooth_seg.csv", index = False)






