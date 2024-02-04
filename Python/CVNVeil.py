import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
    usage='use "python3 %(prog)s --help" for more information',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam_dirs','-bam', help = "Input directory that contains BAM files; if you want to provide multiple directories, separate them by blank space", 
                    nargs= '+')
parser.add_argument('--pre_bam_list','-baml', help = "either bam_dirs or pre_bam_list need to be provided;provide multiple bam files, separate them by blank space", 
                    nargs= '+')
parser.add_argument('--reference','-ref')
parser.add_argument('--reftype','-rt', choices=['hg19','hg38'])
parser.add_argument('--seq_type','-st', choices=['paired-end','single-end'])
parser.add_argument('--output_dir','-o')
parser.add_argument('--prefix','-px',default = "sample")
parser.add_argument('--n_thread','-t',default = 10, type= int)
parser.add_argument('--cell_node','-cn',help = "cell node file; optional, if not given, will not transform cell to node in the end")
parser.add_argument('--change_rate_mode','-crm', choices= ['q','h'], 
    help = "q -> 0.84 quantile; h -> 0.1 hardcode",
    default = 'q' )
parser.add_argument('--ploidy_mode','-pm', choices= ['gl','cl'], 
    help = "global or by cluster",
    default = 'cl' )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
# input_dir = args.input_dir
bam_dirs = args.bam_dirs
pre_bam_list = args.pre_bam_list
reference_file = args.reference
reftype = args.reftype
seq_type = args.seq_type
output_dir = args.output_dir
cell_node = args.cell_node
crm = args.change_rate_mode
pm = args.ploidy_mode
prefix = args.prefix
n_thread = args.n_thread
print(pre_bam_list)
# exit()

import pandas as pd
from collections import Counter
import numpy as np
from sklearn.cluster import AgglomerativeClustering
# import matplotlib.pyplot as plt
import os 
from preprocess import *
from normal_cell_correction import normal_cell_correction

if not os.path.exists(output_dir):
    os.system("mkdir -p "+output_dir)


 
     
def find_local_minima(numbers):
    local_minima = []
    n = len(numbers)

    if n < 3:
        # If the list has less than 3 elements, it cannot have local minima
        return local_minima
    
    idxl = []

    for i in range(1, n - 1):
        if numbers[i] < numbers[i - 1] and numbers[i] < numbers[i + 1]:
            local_minima.append(numbers[i])
            idxl.append(i)

    return idxl


def get_ploidy_from_one_cell_cnv(cnv_oc):
    cand_p =  np.arange(1.15,5.5,0.05)
    rs_list = []
    for p in cand_p:
        rs = abs(cnv_oc*p - np.round(cnv_oc * p)).sum()
        rs_list.append(rs)
    min_idx  = np.argmin(rs_list)
    best_p = cand_p[min_idx]
#     print(min(rs_list),best_p)
    
    idxl  = find_local_minima(rs_list)
    
    lm_p = [ cand_p[i] for i in idxl]
    lm_rs = [rs_list[i] for i in idxl]
    best_p = lm_p[np.argmin(lm_rs)]
#     print(lm_p)
    
    
    # plt.scatter(cand_p,rs_list)
    # plt.vlines(best_p,min(rs_list),max(rs_list),color='red')
    # plt.vlines(lm_p,min(rs_list),max(rs_list),color='green')
    # plt.savefig(output_dir+'/ploidy_'+prefix+'.pdf')
    # # plt.show()
    # plt.close()
    # print("minimum r squre: ",min(lm_rs),'best ploidy: ',best_p)
    
    # print(best_p, min(lm_rs))
    return best_p

def get_ploidy_mt_cell(dat):
    n = len(dat)
    m_list = []
    for i in range(n):
        m = get_ploidy_from_one_cell_cnv(dat[i])
        m_list.append(m)
    print((max(m_list) - min(m_list)))
    print(m_list)
    p_med = np.quantile(m_list,0.5)
    p_mean = get_ploidy_from_one_cell_cnv(dat.mean(0))
    # print((abs(np.array(m_list) - p_med)).mean())
    # if ((max(m_list) - min(m_list))>=1) and ( abs(p_med - p_mean) <= 1.5 ) \
    #     and ((abs(np.array(m_list) - p_med)<=0.2).mean() >=0.95) :
    #     return p_med
    # else:
    #     return p_mean

    return p_mean

def BND(cnv_list):
    return(np.where(cnv_list[1:] - cnv_list[:-1])[0]+0.5)


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
        
def eBND(cnv_list,flanking,thresh):
    dc,_ = nbdiff(cnv_list,flanking)
    bnd = []
    for k,v in dc.items():
        if v>= thresh:
            bnd.append(k)
    return bnd
            
def cluster(y4,t):
    #---------------------------extract fingerprint bins
    print("---------------------------extract fingerprint bins")
    diff = (y4.max(0) -  y4.min(0)).reshape(-1, 1)

    temp = diff.T.copy()
    thresh = np.quantile(diff,0.9)
    pos = temp[0]<thresh
    keep = temp[0]>=thresh
    print('proportion of finger print bins:',(temp[0]>thresh).mean())



    #---------------------------hierachical clustering
    print("---------------------------hierachical clustering")
    data = y4[:,keep]
    # hierarchical_cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
    hierarchical_cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward') #,distance_threshold=t)
    labels = hierarchical_cluster.fit_predict(data)
    print('prediction:',labels)
    # print('ground truth:',c_list)
    # print('ari:',adjusted_rand_score(labels, c_list))
    # plt.scatter(x, y, c=labels)
    # plt.show()

    y4_by_cluster = []
    for i in set(labels):
        y4_by_cluster.append(y4[labels==i,:])
    num_cell_clusters = len(set(labels))
    print("num cell clusters: ", num_cell_clusters)
    print(Counter(labels))
    return np.array(labels)

def get_norm_tm_cell(dat_tnbc,t = 20):
    cl_tnbc = cluster(dat_tnbc.T,t)
    normal_cl_idx = []
    tm_cl_idx = []
    for i in range(len(set(cl_tnbc))):
        p = get_ploidy_from_one_cell_cnv(dat_tnbc.T[cl_tnbc==i].mean(0))
        if 1.7 < p <2.2:
            normal_cl_idx.extend( list(np.where(cl_tnbc==i)[0]))
        else:
            tm_cl_idx.extend( list(np.where(cl_tnbc==i)[0]))
    print("num normal", len(normal_cl_idx))
    print("num tm", len(tm_cl_idx))
    return normal_cl_idx,tm_cl_idx

def get_by_cell_cnv(bnd_list, dat, p):
    last_idx  = 0 
    med_dat = []
    bndl= []
    for i in range(len(bnd_list)):

        cidx = bnd_list[i]
        med = np.quantile(dat[last_idx:cidx],0.5)
        med_dat.extend([med]*(cidx - last_idx))
        bndl.append(cidx-last_idx)
        last_idx = cidx

    med_dat = np.array(med_dat)   
    bndl = np.array(bndl)

    #---------------estimate CNV
    est_cnv = (np.round(p*med_dat)).astype(int)
    est_bnd = BND(est_cnv)
    return est_cnv 

def smooth_block(arr):
    med = np.quantile(arr,0.5)
    med_arr = np.array([[med]*arr.shape[1]]*arr.shape[0])
    x = (arr - med_arr) / 2
    y = med_arr + x
    return y


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




#--------------preprocess

preprocess(bam_dirs,pre_bam_list, reference_file, output_dir, reftype, seq_type ,
            n_thread= n_thread)

#--------------normal cell correction
normal_cell_correction(output_dir)

#----------------para  

pdir= output_dir +'/'
rdmatrix = pdir + "RDmatrix.csv"
genome_cov = pdir +"genome_cov.bed"


# check_cn = glob.glob(pdir +"cell_node*")

# if check_cn:
#     cell_node = check_cn[0]
# else:
#     cell_node = None

print("Using input files:")
print(rdmatrix)
print(genome_cov)

if cell_node is not None:
    print(cell_node)


# outfile = pdir + "CNV.csv"





# -----------------load data
# result = pyreadr.read_r(rdmatrix)
# df=result[None]

df = pd.read_csv(rdmatrix)



# replace df column names
df1 = pd.read_csv(genome_cov, sep = '\t')
df.columns = df1.columns.tolist()[3:]

# get ref bins
ref_bins = df1.apply(lambda x: x.chr + ':'+str(x.start)+'-'+str(x.end) , axis =1).values    

dc_ref = {}
for i in range(df1.shape[0]):
    dc_ref[(df1['chr'][i],  df1['end'][i])] = i



Y = df.values
# Add 1 to every element in the matrix Y
Y_plus_1 = Y + 1

# Calculate the column means of Y_plus_1
column_means = np.mean(Y_plus_1, axis=0)

# Normalize each column of Y_plus_1 by dividing it by its respective column mean
Y_norm = Y_plus_1 / column_means

# Calculate the column std of normal
cell_std = np.std(Y_norm, axis=0)
# std_thresh = np.quantile(cell_std,0.5)
std_thresh = 0.25
print(f"Use standard deviation threshold {std_thresh} to separate normal and cancer cells.")


# plt.hist(cell_std, bins = 50)
# plt.savefig(output_dir+f'/cell_std_{prefix}.pdf')
# plt.close()
all_cells = np.array([x.split('.')[0] for x in df1.columns.tolist()[3:]])
n_cell_idx = np.where(cell_std <= std_thresh)[0]
c_cell_idx = np.where(cell_std > std_thresh)[0]


if len(n_cell_idx)<20:
    n_cell_idx, c_cell_idx = get_norm_tm_cell(Y_norm,t = 50)


Y_norm_n_cell = Y_norm[:,n_cell_idx]
Y_norm_c_cell = Y_norm[:,c_cell_idx]
normal_cells =  all_cells[n_cell_idx]
cancer_cells =  all_cells[c_cell_idx]

print(f"Got {len(n_cell_idx)} normal cells, {len(c_cell_idx)} cancer cells")
print("normal cells:")
print(normal_cells)
print("cencer cells:")
print(cancer_cells)

# ---------------deal with cancer cell
print("Processing cancer cells...")
df4 = pd.DataFrame(Y_norm_c_cell)
df4.columns = cancer_cells



y4 = Y_norm_c_cell.T

    

    
#---------------------extract fingerprint features

diff = (y4.max(0) -  y4.min(0)).reshape(-1, 1)

temp = diff.T.copy()
thresh = np.quantile(diff,0.9)
print("bin variance thresh:",thresh)
pos = temp[0]<thresh
keep = temp[0]>=thresh
print(f"{(temp[0]>=thresh).sum()} out of {len(diff)} bins are extracted as fingerprint")

data = y4[:,keep]
data_norm = data/data.mean(1).reshape(-1,1)


#---------------------- hierachical clustering for ploidy

cell_cl_dist = 30

# hierarchical_cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
print(f"Use threshold cell distance threshold {cell_cl_dist} to perform cell clustering")
hierarchical_cluster = AgglomerativeClustering(n_clusters=None, affinity='euclidean', linkage='ward',distance_threshold=cell_cl_dist)
labels = hierarchical_cluster.fit_predict(data_norm)
print("clustering result:\n",labels)
num_cell_clusters = len(set(labels))
print("num cell clusters: ", num_cell_clusters)

y4_by_cluster = []
for i in set(labels):
    y4_by_cluster.append(y4[labels==i,:])

labels  = np.array(list(labels) + [len(set(labels))] * len(n_cell_idx))
y4_by_cluster.append(Y_norm_n_cell.T)
num_cell_clusters = len(set(labels))
print("(adding normal cell) num cell clusters: ", num_cell_clusters)


#---------------------get global ploidy

# if pm == 'gl':
print("Estimate cancer cell global ploidy")
global_p = get_ploidy_from_one_cell_cnv(np.quantile(y4,0.5,0))

#------------------------------CNV estimation
CNV = []
new_cell_list = []
all_cells = np.array(list(cancer_cells) + list(normal_cells))
final_p_list = np.array([-1.1]*len(labels))
for j in range(len(set(labels))):
# j = 0
    # node = cell_list[list(labels).index(j)]
#     node_bnd = dc_bnd_idx[node]
    n  = y4_by_cluster[j].shape[0]
    # print(labels)
    # print(all_cells)
    
    new_cells = all_cells[labels == j]
    new_cell_list.extend(new_cells)
    print("-------------Process cluster ",j+1,f"({n} cells)")

    # ----------------- get bin consensus
    Y_one_cl = y4_by_cluster[j]
    Y_mean = Y_one_cl.mean(0).reshape(1,-1)
    

    dat = Y_mean[0]

    # ------------- get changing rate
    f = 6
    diffl = nbdiff(dat,f)

    #--------change rate cut-off
    if crm == 'q':
        if n<5:
            t = np.round(np.quantile(diffl,0.84),2)
        else:
            t = np.round(np.quantile(diffl,0.84),2)

    elif crm == 'h':
        t = 0.1
    else:
        print("change_rate_mode can only be 'q' or 'h'")
        exit()
    
    
    '''  t：
    0.84 Q

    or 0.1 hard code ???'''

#     t = 0.1
    print("changing rate cut off ",np.round(t,2))

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




    # ------------------ estimate ploidy

    if j == (len(set(labels))-1):
        ## last one is normal cell
        print(f"Estimate normal cell cluster {j} ploidy")
        # p = get_ploidy_mt_cell(Y_one_cl)
        p = 2.0
        print(p)

    elif n >= 5:
        print(f"Estimate cancer cell cluster {j} ploidy")
        # p = get_ploidy_from_one_cell_cnv(Y_mean,'TMcl'+str(j)+'_'+str(n)+'cells')
        p = get_ploidy_mt_cell(Y_one_cl)
        print(p)


    else:
        print(f"small cancer cell cluster {j} ({n} cells), use global ploidy {global_p}")
        p = global_p

    final_p_list[labels==j] = p



print("The ploidy list")
print(final_p_list)


print("-----------------------------------------------------\n\n")


print("                   Second round clustering           \n\n")


print("-----------------------------------------------------\n\n")

#---------------------- hierachical clustering

cell_cl_dist = 14
# hierarchical_cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
print(f"Use cell distance threshold {cell_cl_dist} to perform cell clustering")
hierarchical_cluster = AgglomerativeClustering(n_clusters=None, affinity='euclidean', linkage='ward',
                                               distance_threshold=cell_cl_dist)
labels = hierarchical_cluster.fit_predict(data_norm)
print("clustering result:\n",labels)
num_cell_clusters = len(set(labels))
print("num cell clusters: ", num_cell_clusters)

y4_by_cluster = []
for i in set(labels):
    y4_by_cluster.append(y4[labels==i,:])

labels  = np.array(list(labels) + [len(set(labels))] * len(n_cell_idx))
y4_by_cluster.append(Y_norm_n_cell.T)
num_cell_clusters = len(set(labels))
print("(adding normal cell) num cell clusters: ", num_cell_clusters)


#---------------------get global ploidy

# if pm == 'gl':
# print("Estimate cancer cell global ploidy")
# global_p = get_ploidy_from_one_cell_cnv(np.quantile(y4,0.5,0))

#------------------------------CNV estimation
CNV = []
new_cell_list = []
all_cells = np.array(list(cancer_cells) + list(normal_cells))
for j in range(len(set(labels))):
# j = 0
    # node = cell_list[list(labels).index(j)]
#     node_bnd = dc_bnd_idx[node]
    Y_one_cl = y4_by_cluster[j]
    n  = y4_by_cluster[j].shape[0]
    # print(labels)
    # print(all_cells)
    
    new_cells = all_cells[labels == j]
    new_cell_list.extend(new_cells)
    print("-------------Process cluster ",j+1,f"({n} cells)")

    # ----------------- get bin consensus
    Y_mean = y4_by_cluster[j].mean(0).reshape(1,-1)

    dat = Y_mean[0]

    # ------------- get changing rate
    f = 6
    diffl = nbdiff(dat,f)

    #--------change rate cut-off
    if crm == 'q':
        if n<5:
            t = np.round(np.quantile(diffl,0.84),2)
        else:
            t = np.round(np.quantile(diffl,0.84),2)

    elif crm == 'h':
        t = 0.1
    else:
        print("change_rate_mode can only be 'q' or 'h'")
        exit()
    
    
    '''  t：
    0.84 Q

    or 0.1 hard code ???'''

#     t = 0.1
    print("changing rate cut off ",np.round(t,2))

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



    for m in range(Y_one_cl.shape[0]):
        p_use = final_p_list[idx_list[m]]
        print("p_use:",p_use)
        est_cnv = get_by_cell_cnv(bnd_list, Y_one_cl_sm[m], p_use)
        # print(Y_one_cl[m].shape,Y_mean[0].shape)
        # est_cnv = get_by_cell_cnv(bnd_list, Y_mean[0], p_use)
        CNV.append(est_cnv)






nbin = len(CNV[0])
# CNV.extend([[2]*nbin]*len(normal_cells))

if cell_node is not None:
    df3 = pd.read_csv(cell_node, sep = '\t', header = None)
    dc_cn = dict(zip(df3[0],df3[1]))
    cells = [dc_cn[cell]+'_'+cell if cell in dc_cn else cell for cell in new_cell_list ]
    # print(transformed_names)
    # print(normal_cells)
    # cells = transformed_names + list(normal_cells)
else:
    # cells = new_cell_list + list(normal_cells)
    cells = new_cell_list

#---------------write to CSV
CNV = np.array(CNV)
# bins = df['bin'].values
bins = df1.apply(lambda x: x.chr+':'+str(x.start)+'-'+str(x.end), axis =1).values
dfcnv = pd.DataFrame(CNV, columns =  bins)
dfcnv['cell'] = cells
dfcnv = dfcnv.set_index('cell')
output_path = output_dir+f'/CNV_{prefix}.csv'
dfcnv.to_csv(output_path)

dfcnv1 = pd.DataFrame(CNV.T, columns =  cells )
dfcnv1['bin'] = bins
dfcnv1 = dfcnv1.set_index('bin')
dfcnv1.to_csv(output_path+'.T')

