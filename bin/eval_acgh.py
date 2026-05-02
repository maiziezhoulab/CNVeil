import pandas as pd
import sys
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy import stats
from tqdm import tqdm
T10_acgh = pd.read_csv("/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/pipeline_test/Real_Data/T10_CN_all.tsv", sep="\t")
genome = T10_acgh[['CHROM', 'START', 'END']]
print(genome.head())
print(T10_acgh.head())
# get cell 2 sector information
# when a clone is found in multple sector, the use the sector with min error 
cell2sec = pd.read_csv("/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/pipeline_test/Real_Data/cell2sector", sep=" ", header=None)
cell2sec =  dict(zip(cell2sec[0], cell2sec[1]))
cn = pd.read_csv(sys.argv[1])

cols = cn.columns.tolist()
if cols[:3]!= ['CHROM','START','END']:
  print(cn.head())
  cells = cols[1:]
  bcol = cols[0]
  bp = cn.iloc[:,0][0]
  if '-' in bp:
    cl = cn[bcol].apply(lambda x: x.split(':')[0])
    sl = cn[bcol].apply(lambda x: int(x.split(':')[1].split('-')[0]))
    el = cn[bcol].apply(lambda x: int(x.split(':')[1].split('-')[1]))
  else:
    cl = cn[bcol].apply(lambda x: x.split(':')[0])
    sl = cn[bcol].apply(lambda x: int(x.split(':')[1]))
    el = cn[bcol].apply(lambda x: int(x.split(':')[2]))

  cn['CHROM'] = cl 
  cn['START'] = sl 
  cn['END'] = el 

  cn = cn[['CHROM','START','END'] + cells]



print(cn.head())
D_err = []
H_err = []
A1_err = []
A2_err = []
newdf = []
for idx, row in genome.iterrows():
  chrom = row['CHROM']
  start = row['START']
  end = row['END']
  # Filter rows in cn that meet the conditions
  temp_df = cn[(cn['CHROM'] == chrom) & (cn['START'] <= start) & (cn['END'] >= end)]
  # Proceed if any rows meet the condition
  if not temp_df.empty:
    temp = temp_df.iloc[0].tolist()  # only 1 possible match
    temp[1] = start  # Update the start value
    temp[2] = end  # Update the end value
    newdf.append(temp)
# Create the new DataFrame
newdf = pd.DataFrame(newdf, columns=cn.columns)
# Create a new dataframe from the filtered dataframe
genome = newdf[['CHROM', 'START', 'END']]
T10_acgh.drop_duplicates(inplace=True)
print(T10_acgh.shape)
acgh = genome.merge(T10_acgh, how="inner")
acgh.dropna(inplace=True)
print(newdf.shape, genome.shape, acgh.shape)
D = []
A1 = []
H = []
A2 = []
D_p = []
A1_p = []
H_p = []
A2_p = []
for cell in tqdm(newdf.columns.to_list()[3:]):
  if cell in cell2sec:
    sec1 = cell2sec[cell]
    #cell_cn = np.log2(0.000001 + newdf[cell]/newdf[cell].mean())
    cell_cn = newdf[cell]
    if sec1.endswith("D"):
      err = mean_squared_error(cell_cn,round(2**(acgh[sec1])*2))
      pearson_coef = stats.pearsonr(cell_cn,round(2**(acgh[sec1])*2))[0]
      print(cell, sec1,cell_cn,round(2**(acgh[sec1])*2) )
      D.append(err)
      D_p.append(pearson_coef)
    elif sec1.endswith("1"):
      err = mean_squared_error(cell_cn, round(2**(acgh[sec1])*3.05))
      A1.append(err)

      pearson_coef = stats.pearsonr(cell_cn,round(2**(acgh[sec1])*3.05))[0]
      A1_p.append(pearson_coef)
    elif sec1.endswith("2"):
      err = mean_squared_error(cell_cn, round(2**(acgh[sec1])*2.85))
      A2.append(err)

      pearson_coef = stats.pearsonr(cell_cn,round(2**(acgh[sec1])*2.85))[0]
      A2_p.append(pearson_coef)
    else:
      err = mean_squared_error(cell_cn, round(2**(acgh[sec1])*1.7))
      H.append(err)

      pearson_coef = stats.pearsonr(cell_cn,round(2**(acgh[sec1])*1.7))[0]
      H_p.append(pearson_coef)
  else:
    print(cell, " not in cell2sec ")
outfile = open(sys.argv[2], "w")

# print("D H A")
de = sum(D)/len(D)
he = sum(H)/len(H)
a1e = sum(A1)/len(A1)
a2e = sum(A2)/len(A2)
print("D error: ", de)
print("H error: ", he)
print("A1 error: ", a1e)
print("A2 error: ", a2e)

outfile.write(f"D_error{de}," + ",".join(str(e) for e in D)+"\n")
outfile.write(f"h_error{he}," + ",".join(str(e) for e in H)+"\n")
outfile.write(f"A1_error{a1e}," + ",".join(str(e) for e in A1)+"\n")
outfile.write(f"A2_error{a2e}," + ",".join(str(e) for e in A2)+"\n")

#-------------print pearson corr
dp = sum(D_p)/len(D_p)
hp = sum(H_p)/len(H_p)
a1p = sum(A1_p)/len(A1_p)
a2p = sum(A2_p)/len(A2_p)
print("D Pearson correlation coefficient : ", dp)
print("H Pearson correlation coefficient : ", hp)
print("A1 Pearson correlation coefficient : ", a1p)
print("A2 Pearson correlation coefficient : ", a2p)

outfile.write(f"D_pearson{dp}," + ",".join(str(e) for e in D_p)+"\n")
outfile.write(f"h_pearson{hp}," + ",".join(str(e) for e in H_p)+"\n")
outfile.write(f"A1_pearson{a1p}," + ",".join(str(e) for e in A1_p)+"\n")
outfile.write(f"A2_pearson{a2p}," + ",".join(str(e) for e in A2_p)+"\n")
