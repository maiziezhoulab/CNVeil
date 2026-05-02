import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
    usage='use "python3 %(prog)s --help" for more information',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--merged_bam_file','-mgbam')
parser.add_argument('--split_bam_dir','-spbam')
parser.add_argument('--reftype','-rt', choices=['hg19','hg38'])
parser.add_argument('--seq_type','-st', choices=['paired-end','single-end'])
parser.add_argument('--data_type','-dtype', choices=['10x','other'])
parser.add_argument('--work_dir','-w')
parser.add_argument('--num_threads','-t', type = int, default= 20)
parser.add_argument('--sample','-sp',default = "sample")
parser.add_argument('--keep_low_ploidy_cell','-kl', action='store_true')
parser.add_argument('--cell_node','-cn',help = "this is designed for the simulation data, not meaningful for real data. cell node file; optional, if not given, will not transform cell to node in the end")
parser.add_argument('--change_rate_mode','-crm', choices= ['q','h'], 
    help = "q -> 0.84 quantile; h -> 0.1 hardcode",
    default = 'q' )
parser.add_argument('--ploidy_mode','-pm', choices= ['gl','cl'], 
    help = "global or by cluster",
    default = 'cl' )

args = parser.parse_args()

merged_bam_file = args.merged_bam_file
split_bam_dir = args.split_bam_dir
reftype = args.reftype
seq_type = args.seq_type
output_dir = args.work_dir + "/Total_CN/"
cell_node = args.cell_node
crm = args.change_rate_mode
pm = args.ploidy_mode
prefix = args.sample
num_threads = args.num_threads
dtype = args.data_type
keep_low_ploidy_cell = args.keep_low_ploidy_cell

import pandas as pd
from collections import Counter, defaultdict
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import os 
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.mixture import GaussianMixture
import pickle
from Perform_QC import perform_qc
from subprocess import Popen
import os



# import matplotlib.pyplot as plt
if not os.path.exists(output_dir):
    os.system("mkdir -p "+output_dir)




def bin_reads_GC(merged_bam, split_bam_folder, refv, out_dir,cell_node_file, num_threads):
    ''' refv : hg19|hg38 '''
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    cmd = f"Rscript {SCRIPT_DIR}/Preprocess/preprocess.R {merged_bam} {split_bam_folder} {refv} {out_dir}/bin_reads/  {num_threads}"
    Popen(cmd, shell= True).wait()

    if cell_node_file is None:
        cmd = f"Rscript {SCRIPT_DIR}/Preprocess/reformat_bin_reads.R  -i {out_dir}/bin_reads/binned-GC  -o {out_dir}/RDmatrix.csv"
        Popen(cmd, shell= True).wait()
    else:
        cmd = f"Rscript {SCRIPT_DIR}/Preprocess/reformat_bin_reads.R  -i {out_dir}/bin_reads/binned-GC  -o {out_dir}/RDmatrix.csv -t {cell_node_file}"
        Popen(cmd, shell= True).wait()




def nbdiff(cnv_list,flanking):
#     dc ={}
    diffl = [0]*flanking
    
    for i in range(flanking, len(cnv_list) - flanking):
        lmed = np.quantile(cnv_list[i-flanking:i],0.5)
        rmed = np.quantile(cnv_list[i:i+flanking],0.5)
        diff = abs(lmed-rmed)
#         dc[i-0.5] = diff 
        diffl.append(diff)
    diffl.extend([0]*flanking)
    return np.array(diffl)
        

            



def smooth_block(arr):
    med = np.quantile(arr,0.5)
    med_arr = np.array([[med]*arr.shape[1]]*arr.shape[0])
    x = (arr - med_arr) / 2
    y = med_arr + x
    return y
    # #--------------- disable smoothing effect
    # return arr


def smooth_cnv(bnd_list, dat):
    last_idx  = 0 
    med_dat = []
    bndl= []
    result = []
    for i in range(len(bnd_list)):

        cidx = bnd_list[i]
        
        x = dat[:,last_idx:cidx]
        y = smooth_block(x)
        if i ==0:
            result = y
        else:
            result = np.column_stack((result,y))
        last_idx = cidx
    return result

def sose(scnp):
    rcnp = np.round(scnp)

    return np.sum(np.square(rcnp-scnp))

def average_seg(boundaries, chrom_matrix):
    # chrom_matrix = chrom_matrix.T
    for i in range(len(chrom_matrix)):
        for j in range(len(boundaries)):
            if j == 0:
                chrom_matrix[i][0:boundaries[j]] = np.median(chrom_matrix[i][0:boundaries[j]])
            else:
                chrom_matrix[i][boundaries[j-1]:boundaries[j]] = np.median(chrom_matrix[i][boundaries[j-1]:boundaries[j]])
    return chrom_matrix

def average_seg_multi_cell(boundaries, chrom_matrix):
    # Iterate through boundary segments
    for j in range(len(boundaries)):
        if j == 0:
            # Get the segment for all cells
            segment = chrom_matrix[:, 0:boundaries[j]]
            # Calculate the median across all cells
            segment_median = np.median(segment)
            # Assign the median back to all cells in the segment
            chrom_matrix[:, 0:boundaries[j]] = segment_median
        else:
            # Get the segment for all cells
            segment = chrom_matrix[:, boundaries[j-1]:boundaries[j]]
            # Calculate the median across all cells
            segment_median = np.median(segment)
            # Assign the median back to all cells in the segment
            chrom_matrix[:, boundaries[j-1]:boundaries[j]] = segment_median
    return chrom_matrix


def estimate_candidate_ploidy(Y):
    # Y = Y.T
    Y_nor = []
    min_ploidy = 1.0
    max_ploidy = 5.5
    X = np.arange(min_ploidy, max_ploidy, 0.05)
    local_min_p_list = []
    local_min_sos_list = []
    p_list = []
    r_list = []
    smaller_cand_p = []
    bigger_cand_p = []
    
    for cell in Y:
        RCNP = cell / np.mean(cell)
        sose_list = []
        for i in X:
            sose_list.append(sose(RCNP * i))
        

        # Collect local minima and corresponding X (ploidy)
        local_min_indices = (np.diff(np.sign(np.diff(sose_list))) > 0).nonzero()[0] + 1
        local_min_sos = [sose_list[idx] for idx in local_min_indices]
        local_min_p = [X[idx] for idx in local_min_indices]

        
        ploidy = round(X[sose_list.index(min(sose_list))],2)
        if ploidy==min(X):
            ploidy = round(min(local_min_p),2)
        # print("get ploidy: ", ploidy)
        p_list.append(ploidy)
        SCNP = list(np.round(RCNP * ploidy))
        Y_nor.append(SCNP)

        sorted_min_sos = sorted(local_min_sos)
        r = sorted_min_sos[1]/sorted_min_sos[0]

        top2_indices = np.argsort(local_min_sos)[:2]
        top2_best_p = sorted([local_min_p[top2_indices[0]],local_min_p[top2_indices[1]]])
        smaller_cand_p.append(top2_best_p[0])
        bigger_cand_p.append(top2_best_p[1])

        local_min_sos_list.extend(local_min_sos)
        local_min_p_list.extend(local_min_p)
        r_list.append(r)
        



    return Y_nor, local_min_p_list, local_min_sos_list,p_list, np.array(r_list), smaller_cand_p, bigger_cand_p


def call_CN(Y,ploidy):
    Y_nor = []
    for i,cell in enumerate(Y):
        RCNP = cell / np.mean(cell)
        if type(ploidy) is list:
            use_p = ploidy[i]
        else:
            use_p = ploidy
        SCNP = list(np.round(RCNP * use_p))
        Y_nor.append(SCNP)
    return np.array(Y_nor)








def filter_low_cov_bins(Y,min_cov ):
    # Compute normalization factors
    bin_filter = np.apply_along_axis(lambda x: x.mean()>min_cov , axis=1, arr=Y)
    # filtered_rows = np.apply_along_axis(lambda x: np.sum(x == 0) < 0.1 * len(x), axis=1, arr=Y)
    n = len(bin_filter) - bin_filter.sum()
    print(f"{n} bins out of {Y.shape[0]} bins will be removed due to low coverage")
    print(f"{Y.shape[0] - n} bins will be left")
    return bin_filter



def preprocess(output_dir, reftype, dtype , keep_low_ploidy_cell):

    #----------------para  

    pdir= output_dir +'/'
    rc_file = 'RDmatrix_qc.csv'
    rdmatrix = pdir + rc_file


    print("Using input files:")
    print(os.path.join(pdir,'RDmatrix.csv') )

    if cell_node is not None:
        print("Use cell_node converter: ", cell_node)
        
    # ----------step1 QC
    # if not os.path.exists(rdmatrix):
    if dtype == '10x':
        perform_qc(pdir+"/RDmatrix.csv", reftype, rdmatrix, False, keep_low_ploidy_cell)
    else:
        perform_qc(pdir+"/RDmatrix.csv", reftype, rdmatrix, True, keep_low_ploidy_cell)

    # load RC matrix
    df = pd.read_csv(rdmatrix)
    rc_mean_by_bin = df.iloc[:,1:].values.mean(axis = 1)


    # -------------step2 normalize RD matrix
    bins = df.iloc[:,0]

    # each column is cell, each row is bin
    Y = df.iloc[:,1:].values

    # Add 1 to every element in the matrix Y
    Y_plus_1 = Y + 1

    # Calculate the column means of Y_plus_1
    column_means = np.mean(Y_plus_1, axis=0)

    # Normalize each column of Y_plus_1 by dividing it by its respective column mean
    Y_norm = Y_plus_1 / column_means

    all_cells = np.array([x.split('.')[0] for x in df.columns.tolist()[1:]])
    # -------------step3  estimate_candidate_ploidy
    # Y_nor, local_min_p_list, local_min_sos_list,p_list, r_list,  smaller_cand_p, bigger_cand_p = estimate_candidate_ploidy(Y_norm.T)
    return Y_norm, all_cells, bins



def identify_normal_cancer_cells(Y_norm, all_cells):
    # -------------step4 Identify normal cells
    # 
    # # Calculate the column std of normal
    cell_std = np.std(Y_norm, axis=0)

    # plot hist
    plt.hist(cell_std, bins = 50)
    plt.savefig(output_dir+"/cell_std.pdf")
    plt.close()

    plt.hist(cell_std[cell_std<1], bins = 50)
    plt.savefig(output_dir+"/cell_std_less_than1.pdf")
    plt.close()

    # std_thresh = np.quantile(cell_std,0.5)
    std_thresh = 0.25
    print(f"Use standard deviation threshold {std_thresh} to separate normal and cancer cells.")


    n_cell_idx = np.where(cell_std <= std_thresh)[0]
    c_cell_idx = np.where(cell_std > std_thresh)[0]



    Y_norm_n_cell = Y_norm[:,n_cell_idx]
    Y_norm_c_cell = Y_norm[:,c_cell_idx]
    normal_cells =  all_cells[n_cell_idx]
    cancer_cells =  all_cells[c_cell_idx]

    print(f"Got {len(n_cell_idx)} normal cells, {len(c_cell_idx)} cancer cells")
    print("normal cells:")
    print(normal_cells)
    print("cencer cells:")
    print(cancer_cells)
    return Y_norm_n_cell, Y_norm_c_cell, normal_cells, cancer_cells, n_cell_idx, c_cell_idx


def cluster_prepare(Y_norm_c_cell,r1,r2):
    # ---------------deal with cancer cell
    # now each row is cell
    y4 = Y_norm_c_cell.T

    #---------------------extract fingerprint features

    diff = (y4.max(0) -  y4.min(0)).reshape(-1, 1)
    temp = diff.T.copy()
    thresh = np.quantile(diff,0.9)
    print("Use 90% quantile of bin variance as thresh for fingerprint bins\nThreshold=",thresh)
    pos = temp[0]<thresh
    keep = temp[0]>=thresh
    print(f"{(temp[0]>=thresh).sum()} out of {len(diff)} bins are extracted as fingerprint")
    data = y4[:,keep]
    #normalize data again using only fingerprint bins
    data_norm = data/data.mean(1).reshape(-1,1)

    #---------------------- hierachical clustering for ploidy
    # Create the linkage matrix for the dendrogram
    linkage_matrix = linkage(data_norm, method='ward')
    # Plot the dendrogram
    plt.figure(figsize=(10, 7))
    dendrogram(linkage_matrix)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Samples')
    plt.ylabel('Distance')
    plt.savefig(output_dir+"/dendrogram.pdf")
    plt.close()

    # Get the top height of the dendrogram
    top_height = round(np.max(linkage_matrix[:, 2]))
    print(f"The top height of the dendrogram is: {top_height}")
    cell_cl_dist1 = round(r1 * top_height)
    cell_cl_dist2 = round(r2 * top_height)
    print(f"First round clustering dist threshold: {r1} * {top_height}(top height) = {cell_cl_dist1}")
    print(f"Second round clustering dist threshold: {r2} * {top_height}(top height) = {cell_cl_dist2}")
    with open(output_dir + "/cell_cluster_dist",'w' ) as f:
        f.write(f"Dendrogram top height: {top_height}\n")
        f.write(f"First round clustering distance threshold({r1} * top height): {cell_cl_dist1}\n")
        f.write(f"Second round clustering distance threshold({r2} * top height): {cell_cl_dist2}\n")

    return data_norm, top_height, cell_cl_dist1, cell_cl_dist2


def first_round_cluster(Y_norm_c_cell, data_norm,  n_cell_idx,cell_cl_dist1):

    print("-----------------------------------------------------\n\n")


    print("                   First round clustering           \n\n")


    print("-----------------------------------------------------\n\n")
    y4 = Y_norm_c_cell.T
    print(f"Use cell distance threshold {cell_cl_dist1} to perform first round cell clustering")
    hierarchical_cluster = AgglomerativeClustering(n_clusters=None, metric='euclidean', linkage='ward',distance_threshold=cell_cl_dist1)
    labels = hierarchical_cluster.fit_predict(data_norm)
    # print("clustering result:\n",labels)
    num_cell_clusters = len(set(labels))
    print("num cell clusters (cancer cells only): ", num_cell_clusters)
    print(Counter(labels))


    #  each element is a cluster; within cluster, each row is cell
    y4_by_cluster = []
    cells_by_cluster = []
    for i in set(labels):
        y4_by_cluster.append(y4[labels==i,:])
        cells_by_cluster.append(cancer_cells[labels==i])
    # print("len(cells_by_cluster)",len(cells_by_cluster))
    first_round_dict_label2cell = defaultdict(list)
    for i in range(len(labels)):
        label = labels[i]
        cell = cancer_cells[i]
        first_round_dict_label2cell[label].append(cell)

    # Save the dictionary to a pickle file
    with open(output_dir+"/ploidy_cluster_result.pkl" , "wb") as f:  # Use 'wb' mode for writing in binary
        pickle.dump(first_round_dict_label2cell, f)
    # print(first_round_dict_label2cell)
    # exit()

    labels  = np.array(list(labels) + [len(set(labels))] * len(n_cell_idx))



    y4_by_cluster.append(Y_norm_n_cell.T)
    cells_by_cluster.append(normal_cells)

    num_cell_clusters = len(set(labels))
    print("(adding normal cell) num cell clusters: ", num_cell_clusters)
    return first_round_dict_label2cell


def second_round_cluster(Y_norm_c_cell,data_norm,first_round_dict_label2cell,cell_cl_dist2):
    print("-----------------------------------------------------\n\n")


    print("                   Second round clustering           \n\n")


    print("-----------------------------------------------------\n\n")

    #---------------------- hierachical clustering
    y4 = Y_norm_c_cell.T
    # hierarchical_cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
    print(f"Use cell distance threshold {cell_cl_dist2} to perform cell clustering")
    hierarchical_cluster = AgglomerativeClustering(n_clusters=None, metric='euclidean', linkage='ward',
                                                distance_threshold=cell_cl_dist2)
    labels = hierarchical_cluster.fit_predict(data_norm)
    print("clustering result:\n",labels)
    num_cell_clusters = len(set(labels))
    print("num cell clusters: ", num_cell_clusters)
    second_round_dict_cell2label = dict(zip(cancer_cells,labels))
    # print(second_round_dict_cell2label)
    # connect first and second round clustering
    first2second_dict = {}
    for label, cells in first_round_dict_label2cell.items():
        # print([second_round_dict_cell2label[cell] for cell in cells])
        first2second_dict[label] = list(set([second_round_dict_cell2label[cell] for cell in cells]))

    print(first2second_dict)
    # exit()
    y4_by_cluster = []
    for i in set(labels):
        y4_by_cluster.append(y4[labels==i,:])

    labels  = np.array(list(labels) + [len(set(labels))] * len(n_cell_idx))

    for cell in normal_cells:
        second_round_dict_cell2label[cell] = len(set(labels))

    second_round_dict_label2cell = defaultdict(list)
    for cell, id in second_round_dict_cell2label.items():
        second_round_dict_label2cell[id].append(cell)
    # print(cancer_cells)
    # print(normal_cells)
    # print(labels)
    # second_round_dict_label2cell = dict(zip(labels, list(cancer_cells)+list(normal_cells)))
    with open(output_dir+"/fine_cluster_result.pkl" , "wb") as f:  # Use 'wb' mode for writing in binary
        pickle.dump(second_round_dict_label2cell, f)
    # exit()
    y4_by_cluster.append(Y_norm_n_cell.T)
    num_cell_clusters = len(set(labels))
    print("(adding normal cell) num cell clusters: ", num_cell_clusters)
    return first2second_dict, labels, y4_by_cluster

def double_clustering(Y_norm_c_cell,  r1, r2 ):
    data_norm, top_height, cell_cl_dist1, cell_cl_dist2 = \
        cluster_prepare(Y_norm_c_cell, r1, r2)

    first_round_dict_label2cell \
        = first_round_cluster(Y_norm_c_cell, data_norm,  n_cell_idx, cell_cl_dist1)

    first2second_dict, second_cluster_labels, RDmatrix_by_second_cluster \
        = second_round_cluster(Y_norm_c_cell, data_norm, first_round_dict_label2cell,cell_cl_dist2)
    return RDmatrix_by_second_cluster, second_cluster_labels ,first2second_dict

def estimate_segmentation(y4_by_cluster, labels):
    print("-----------------------------------------------------\n\n")


    print("                   Estimate segmentation           \n\n")


    print("-----------------------------------------------------\n\n")
    #------------------------------CNV estimation
    CNV = []
    new_cell_list = []
    new_cell_list1 = []
    all_cells = np.array(list(cancer_cells) + list(normal_cells))
    print("all_cells",all_cells[:10])
    # flog = open(output_dir+"/ploidy_by_cell_by_cluster.txt",'w')
    local_min_p_list, local_min_sos_list = [],[]
    avg_seg_list =[]
    avg_seg_list_multi =[]
    p_list_by_cluster = []
    r_list_by_cluster = []
    bnd_all_cells = []
    for j in range(len(set(labels))):

        Y_one_cl = y4_by_cluster[j]

        Y_nor, local_min_p_list, local_min_sos_list,p_list, r_list,  smaller_cand_p, bigger_cand_p = estimate_candidate_ploidy(Y_one_cl)
        p_list_by_cluster.append(p_list)
        r_list_by_cluster.append(r_list)

        # number of cells
        n  = Y_one_cl.shape[0]    
        new_cells = all_cells[labels == j]
        new_cell_list.extend(list(new_cells))
        new_cell_list1.append(list(new_cells))
        # print("-------------Process cluster ",j+1,f"({n} cells)")
        # print(np.median(p_list))

        # ----------------- get bin consensus
        Y_mean = Y_one_cl.mean(0).reshape(1,-1)

        dat = Y_mean[0]

        # ------------- get changing rate
        f = 6
        diffl = nbdiff(dat,f)

        #--------change rate cut-off
        if crm == 'q':
            
            t = np.round(np.quantile(diffl,0.84),2)
            # print(f"Changing rate threshold: 84% quantile of global changing rate {t}")

        elif crm == 'h':
            # print("Changing rate threshold: hard code 0.1")
            t = 0.1
        else:
            print("change_rate_mode can only be 'q' or 'h'")
            exit()
        
        # print("Use changing rate threshold: ",np.round(t,2))

        #---------------get segmentation
        apx_ids = np.where(diffl>= t)[0]

        cl_list = [[apx_ids[0]]]

        for i in range(1,len(apx_ids)):
            cid = apx_ids[i]
            lid = cl_list[-1][-1]
            if cid - lid <=1:
                cl_list[-1].append(cid)
            else:
                cl_list.append([cid])
        clrg_list = [ (min(cl), max(cl)) for cl in cl_list  if len(cl)>=4]

        #---------transform bnd
        a2 = []
        for i in range(len(clrg_list)):
            cl = clrg_list[i]

            a2.append(np.round(cl[0]/2+cl[1]/2)-0.5)
        a2 = np.array(a2)


        bnd_list = list(a2.astype(int)+1)+[len(dat)]

        # ------------------ estimate CNV

        ############ --------------------- do it by cell
        # ---------------get segmentation consensus

        idx_list = np.where(labels == j)[0]
        # Y_one_cl.mean(0)


        Y_one_cl_sm = smooth_cnv(bnd_list, Y_one_cl)

        #---------- test only ----------
        Y_one_cl_avg = average_seg(bnd_list, Y_one_cl)
        Y_one_cl_avg_multi = average_seg_multi_cell(bnd_list, Y_one_cl.copy())
        bnd_all_cells.append((bnd_list,n))

        # est_cnv, local_min_p, local_min_sos,p_list = call_CN(Y_one_cl_avg)
        # flog.write(f"Cluster {j+1} ploidy: "+' '.join([str(p) for p in p_list])+"\n")

        avg_seg_list.append(Y_one_cl_avg)
        avg_seg_list_multi.append(Y_one_cl_avg_multi)
        

    # Save list to pickle file
    with open(output_dir+ "/bnd_all_cells.pkl", "wb") as fb:
        pickle.dump(bnd_all_cells, fb)
    print("write BNDs to "+output_dir+ "/bnd_all_cells.pkl")
    
    return new_cell_list, new_cell_list1, avg_seg_list, avg_seg_list_multi, p_list_by_cluster, r_list_by_cluster



import matplotlib.pyplot as plt
import numpy as np
import math

def plot_ploidy_histograms(cluster_ploidies, outfile="ploidy_hist.png"):
    """
    cluster_ploidies: list of lists
        Each sublist contains ploidy values for one cluster.

    outfile: str
        Output PNG file name.
    """

    num_clusters = len(cluster_ploidies)
    cols = 3  # number of subplots per row
    rows = math.ceil(num_clusters / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes = axes.flatten()  # flatten grid

    for idx, ploidy_list in enumerate(cluster_ploidies):
        ax = axes[idx]

        data = np.array(ploidy_list)

        # statistics
        mean_val = np.mean(data)
        median_val = np.median(data)
        min_val = np.min(data)
        max_val = np.max(data)
        std_val = np.std(data)
        size = len(data)

        # histogram
        ax.hist(data, bins=20)

        # vertical lines for median and spread
        ax.axvline(median_val, color="red", linestyle="--", linewidth=1.5)
        ax.axvline(median_val - std_val, color="green", linestyle=":", linewidth=1.2)
        ax.axvline(median_val + std_val, color="green", linestyle=":", linewidth=1.2)

        # title
        ax.set_title(
            f"Cluster {idx} (n={size})\n"
            f"median={median_val:.2f}, mean={mean_val:.2f}, std={std_val:.2f}"
        )

        # annotation INSIDE the figure
        ax.text(
            0.95, 0.95,
            f"min={min_val:.2f}\nmax={max_val:.2f}",
            transform=ax.transAxes,
            ha='right', va='top',
            fontsize=10,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='none')
        )

    # remove unused axes
    for j in range(num_clusters, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def clamp_list(values, low, high):
    return [max(low, min(v, high)) for v in values]


def estimate_ploidy(y4_by_cluster,first2second_dict,new_cell_list1, avg_seg_list, avg_seg_list_multi, p_list_by_cluster, r_list_by_cluster ):
    print("-----------------------------------------------------\n\n")


    print("                    Estimate ploidy                \n\n")


    print("-----------------------------------------------------\n\n")

    f = open(output_dir+"/mean_ploidy",'w' ) 
    all_cnv = []
    p_list_avg_seg =[]
    p_list_norm_all = []
    final_cell_list = []
    cluster_ploidies = []
    for label, idxs in first2second_dict.items():
        print("-----------process cluster",label)
        avg_seg_one_cluster = []
        avg_seg_one_cluster_multi = []
        cell_name_one_cluster = []
        RD_one_cluster = []
        p_list_one_cluster = []
        r_list_one_cluster = []
        for idx in idxs:
            final_cell_list.extend(new_cell_list1[idx])
            avg_seg_one_cluster.extend(avg_seg_list[idx])
            avg_seg_one_cluster_multi.extend(avg_seg_list_multi[idx])
            RD_one_cluster.extend(y4_by_cluster[idx])
            cell_name_one_cluster.extend(new_cell_list1[idx])
            p_list_one_cluster.extend(p_list_by_cluster[idx])
            r_list_one_cluster.extend(r_list_by_cluster[idx])
        Y_nor, local_min_p_list, local_min_sos_list,p_list, r_list,  smaller_cand_p, bigger_cand_p = estimate_candidate_ploidy(avg_seg_one_cluster)
        Y_nor_norm, local_min_p_list_norm, local_min_sos_list_norm,p_list_norm, r_list_norm,  smaller_cand_p_norm, bigger_cand_p_norm = estimate_candidate_ploidy(RD_one_cluster)
        p_list_norm_all.extend(p_list_norm)
        p_list_avg_seg.extend(p_list)
        cluster_ploidies.append(p_list)
        



        
        p_list_std = round(np.std(p_list),4)
        median_p = np.median(p_list)
        if p_list_std > 0.45:
            best_p = median_p
        else:
            low_b = median_p - p_list_std
            high_b = median_p + p_list_std
            best_p = clamp_list(p_list, low_b, high_b)
            



        if type(best_p) is list:
            summary = f"cluster {label}, {len(avg_seg_one_cluster)} cells, ploidy std {p_list_std}, use each cell's ploidy (clamped:med-std,med+std)({round(min(best_p),2)}-{round(max(best_p),2)})\n"
        else:
            summary = f"cluster {label}, {len(avg_seg_one_cluster)} cells, ploidy std {p_list_std}, use median ploidy {round(best_p,4)}\n"
        print(summary)
        f.write(summary)
        cnv_one_cluster = call_CN( np.array(avg_seg_one_cluster)*0.5 + np.array(avg_seg_one_cluster_multi) *0.5 , best_p )
        all_cnv.append(cnv_one_cluster)

    final_cell_list.extend(list(normal_cells))
    p_list_norm_all = np.array(p_list_norm_all)
    p_list_avg_seg = np.array(p_list_avg_seg)
    np.save(output_dir+"/p_list_RDmatrix.npy"%(label), p_list_norm_all)
    np.save(output_dir+"/p_list_avg_seg.npy"%(label),p_list_avg_seg )

    plot_ploidy_histograms(cluster_ploidies, output_dir+"/cluster_ploidy_hist.png")

    avg_p_avg_seg = round(p_list_avg_seg.mean(),2)
    std_p_avg_seg = round(p_list_avg_seg.std(),2)
    avg_p_rd = round(p_list_norm_all.mean(),2)
    std_p_rd = round(p_list_norm_all.std(),2)

    with open(output_dir+"/avg_ploidy_raw.txt",'w' ) as fr:
        fr.write(f"from RDmatrix:\n avg ploidy {avg_p_rd} std {std_p_rd}\n")
        fr.write(f"from avg seg:\n avg ploidy {avg_p_avg_seg} std {std_p_avg_seg}\n")

    CNV = np.vstack(all_cnv)
    print(CNV.shape)
    print(CNV.mean())
    # print(CNV[:, extra_ids].mean())
    # print(CNV[:, ~np.array(extra_ids)].mean())
    print("Overall ploidy: "+str(np.round(CNV.mean(),3))+"\n")
    f.write("Overall ploidy: "+str(np.round(CNV.mean(),3))+"\n")
    # f.close()

    f.close()

    return CNV, final_cell_list



def write_cnv(final_cell_list, CNV, dtype ):
    if cell_node is not None:
        df3 = pd.read_csv(cell_node, sep = '\t', header = None)
        dc_cn = dict(zip(df3[0],df3[1]))
        cells = [dc_cn[cell]+'_'+cell if cell in dc_cn else cell for cell in final_cell_list ]

    else:
        cells = final_cell_list


    output_path = output_dir+f'/CNV_{prefix}.csv'
    normal_cnv = np.round(Y_norm_n_cell * 2 )
    all_cnv = np.hstack((CNV.T, normal_cnv))

    if dtype == '10x':
        if ((normal_cnv.shape[1])/(all_cnv.shape[1]))>0.1:
            bin_filter1 = filter_low_cov_bins(normal_cnv, min_cov = 1.5)
            bin_filter2 = filter_low_cov_bins(CNV.T, min_cov = 0.5)
            bin_filter = bin_filter1 & bin_filter2
        else:
            bin_filter = filter_low_cov_bins(all_cnv, min_cov = 0.5)



    dfcnv1 = pd.DataFrame(all_cnv, columns = cells)
    dfcnv1['bin'] = bins
    if dtype == '10x':
        dfcnv1 = dfcnv1.loc[bin_filter,:]
    dfcnv1 = dfcnv1.set_index('bin')
    dfcnv1.to_csv(output_path+'.T')




    #---------------write to CSV

    dfcnv = pd.DataFrame(all_cnv.T, columns= bins)
    dfcnv['cell'] = list(cells) 
    dfcnv = dfcnv.set_index('cell')
    if dtype == '10x':
        dfcnv = dfcnv.loc[:,bin_filter]
    output_path = output_dir+f'/CNV_{prefix}.csv'
    dfcnv.to_csv(output_path)

def infer_cnv(RDmatrix_by_second_cluster, second_cluster_labels ,first2second_dict, dtype):
    new_cell_list, new_cell_list1, avg_seg_list, \
        avg_seg_list_multi, p_list_by_cluster, r_list_by_cluster =\
    estimate_segmentation(RDmatrix_by_second_cluster, second_cluster_labels )

    CNV, final_cell_list = \
        estimate_ploidy(RDmatrix_by_second_cluster ,first2second_dict,\
                        new_cell_list1, avg_seg_list, avg_seg_list_multi,\
                              p_list_by_cluster, r_list_by_cluster)

    write_cnv(final_cell_list, CNV, dtype)


if __name__ == '__main__':

    #----------------------------------  Parameters

    r1 = 0.6
    r2 = 0.3

    # #------temp for well dr seq
    # r1 = 0.3
    # r2 = 0.2

    # ---------------- bin reads and GC correction

    bin_reads_GC(merged_bam_file, split_bam_dir, reftype, output_dir, cell_node, num_threads)


    #----------------------------------  Load Data and preprocess

    Y_norm, all_cells, bins = preprocess(output_dir, reftype, dtype , keep_low_ploidy_cell)


    #---------------------------------- Identify normal and cancer cells

    Y_norm_n_cell, Y_norm_c_cell, normal_cells, \
        cancer_cells, n_cell_idx, c_cell_idx = \
            identify_normal_cancer_cells(Y_norm, all_cells)

    #---------------------------------- 2 round clustering 

    RDmatrix_by_second_cluster, second_cluster_labels ,first2second_dict = \
        double_clustering(Y_norm_c_cell,  r1, r2 )


    #---------------------------------- Infer CNV 

    infer_cnv(RDmatrix_by_second_cluster, second_cluster_labels ,first2second_dict, dtype)

