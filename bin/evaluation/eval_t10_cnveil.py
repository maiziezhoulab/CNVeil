import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_dir','-o')

args = parser.parse_args()
comp_csv = args.input_path
out_dir = args.output_dir

import pandas as pd 
import numpy as np
from collections import defaultdict
from sklearn.metrics import mean_squared_error
from scipy import stats
import pickle
import os
# comp_csv = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/Use_AF_matrix/T10/CNV_T10.csv.T"
gold_tsv = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/Data/Real_Data/T10/T10_gold_CNV.tsv"
# out_dir = "CNVeil_T10_eval/"
bin_size = 500e3
p_list = np.array([3.05,2.85, 1.7, 2])
subclone_list = ['A1','A2','H','D']
df_gold = pd.read_csv(gold_tsv, sep = '\t')
df_comp = pd.read_csv(comp_csv)
df_comp['CHROM'] =df_comp['bin'].apply(lambda x: x.split(':')[0])
df_comp['start_i'] =df_comp['bin'].apply(lambda x: int(x.split(':')[1].split('-')[0])//bin_size)
# print(df_comp.columns)
cells = df_comp.columns[1:-2]
mean_ps = df_comp[cells].mean().tolist()
# print(df_comp.head())
# print(df_comp.shape)
# print(df_gold.shape)
df_merge = df_comp.merge(df_gold, on = ['CHROM','start_i'], how='inner')

# print(df_merge)
cell2sector = {}
sector2cell = defaultdict(list)
dc_eval_err = defaultdict(list)
dc_eval_pearson = defaultdict(list)
with open("/data/maiziezhou_lab/CanLuo/Single_Cell_Project/T10_EVALUATION/Data/cell2sector_consensus.pkl", "rb") as f:
    dc_concensus = pickle.load(f)
# print(len(cells),len(mean_ps))
for i in range(len(cells)):
	cell = cells[i]
	if cell in dc_concensus:
		# print(dc_concensus[cell])
		mean_p =  mean_ps[i]
		# print(mean_p)
		j = np.argmin(abs(p_list - mean_p))
		# print(mean_p,j, p_list[j])
		cell2sector[cells[i]] = subclone_list[j]
		sector2cell[subclone_list[j]].append(cells[i])
		cell_cn = df_merge[cell]
		subclone = dc_concensus[cell]
		gold_cn = df_merge[subclone]
		err = mean_squared_error(cell_cn, gold_cn)
		pearson_coef = stats.pearsonr(cell_cn,gold_cn)[0]
		dc_eval_err[subclone].append(err)
		dc_eval_pearson[subclone].append(pearson_coef)

if not os.path.exists(out_dir):
	os.system("mkdir -p " + out_dir)
	
with open(f"{out_dir}/cell2sector.pkl", "wb") as f:
	pickle.dump(cell2sector, f)
# print(dc_eval_err)
# print(dc_eval_pearson)



def write_dict(dc, out_dir, prefix, metric):
	# Save to a pickle file
	with open(f"{out_dir}/{prefix}.pkl", "wb") as f:
		pickle.dump(dc, f)
	with open(f"{out_dir}/{prefix}_summary.txt", 'w') as f:
		for k,v in dc.items():
			mean_val = np.mean(v)
			std = np.std(v)
			line = f"{k}({len(v)} cells):mean {metric} value {round(mean_val,4)}, std {round(std,4)}\n"
			print(line[:-1])
			f.write(line)

write_dict(dc_eval_err, out_dir, "MSE", "MSE")
write_dict(dc_eval_pearson, out_dir, "pearson_coef", "pearson_coef")