
import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
	usage='use "python3 %(prog)s --help" for more information',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i')
parser.add_argument('--output_path','-o')
# parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_path = args.output_path
# n_thread = args.n_thread

import pandas as pd
cell2sec = pd.read_csv("/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/pipeline_test/Real_Data/cell2sector", sep=" ", header=None)
cell2sec.columns = ['cell','subclone']
cell2sec.subclone = cell2sec['subclone'].apply(lambda x : x.split('-')[1])

# dc =  dict(zip(cell2sec[0], cell2sec[1]))
# print(dc )

cell_list_ordered = cell2sec['cell'].tolist()


df = pd.read_csv(input_path)
cell_list_given = set(df.columns[1:].tolist())

cell_list_ordered = [ x for x in cell_list_ordered if x in cell_list_given]

df1 = df[['bin'] + cell_list_ordered]

df1.to_csv(output_path,index = False)



