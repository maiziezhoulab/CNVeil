import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
    usage='use "python3 %(prog)s --help" for more information',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input_path','-i', help = "a csv file, columns = tool,ploidy,eval_dir")
parser.add_argument('--output_dir','-o')
# parser.add_argument('--n_thread','-t', type = int, default = 22 )
# parser.add_argument('--delete_temp_file','-d', action='store_true')
args = parser.parse_args()
input_path = args.input_path
output_dir = args.output_dir
# n_thread = args.n_thread


import os
import pandas as pd
import numpy as np


def collect_data(p, tool_list, eval_dir_list, data_type):
    # tool_list = [
    #              "ginkgo",
    #              "hmmcopy",
    #              "scope",
    #              "secnv"]
    if data_type == 'CNV':
        
        file_list = [ eval_dir_list[i]+'/stats_by_cell_2000000_matchcn.csv' for i in range(len(tool_list))]
    else:
        file_list = [ eval_dir_list[i]+'/stats_by_cell_2000000.csv' for i in range(len(tool_list))]
    precision_list =[]
    recall_list =[]
    f1_list = []
    name_list = []
    for i in range(len(tool_list)):
        tool = tool_list[i]
        file = file_list[i]
        if os.path.exists(file):
            df = pd.read_csv(file)
            precision_list.extend(df['precision'].tolist())
            recall_list.extend(df['recall'].tolist())
            f1 = df['precision'] * df['recall'] * 2 /(df['precision'] + df['recall'])
            f1_list.extend(f1)
            name_list.extend([tool]*df.shape[0])
    dc_recall = {'tool':name_list, 'value': recall_list, }
    dc_precision = {'tool':name_list, 'value': precision_list}
    dc_f1 = {'tool':name_list, 'value': f1_list}
    
    df_recall = create_ordered_df(dc_recall, tool_list)
    df_recall['metric'] = 'Recall'
#     print((df_recall['recall']<0).sum())
    
    df_precision = create_ordered_df(dc_precision, tool_list)
    df_precision['metric'] = 'Precision'
#     print((df_precision['precision']<0).sum())

    df_f1 = create_ordered_df(dc_f1, tool_list)
    df_f1['metric'] = 'F1'
    
    df = pd.concat([df_recall, df_precision,df_f1]).reset_index(drop=True)
    df['Ploidy'] = p
    return df_recall, df_precision,df_f1,df

def create_ordered_df(dc, group_order):
    
    # create sample data
    df = pd.DataFrame(dc)

    # specify the desired order of the groups
    # group_order = [
    #              "ginkgo",
    #              "hmmcopy",
    #              "scope",
    #              "secnv"]

    # create a categorical variable for the group column
    df['Tool'] = pd.Categorical(df['tool'], categories=group_order)
    
    return df


def met_stats(df,tool,p,met):
    val = df[((df['Tool'] == tool) & (df['Ploidy']==p) & (df['metric']==met))]['value']
    info = [p,met,tool,len(val),
    val.mean(), 
    np.std(val), 
    val.max(),
    val.min(),
    val.max() - val.min()]

    for i in range(3,len(info)):
        info[i] = round(info[i],3)
    return info


def stats_one_tool(df, tool,p):
    num_cells = int(((df['Tool'] == tool) & (df['Ploidy']==p)).sum()/3)
    infos = [ met_stats(df,tool,p,met)  for met in ['Recall','Precision','F1']]
    # if num_cells:
    return infos

def collect(data_type, tool_list, p_list, df_config, out_dir):
    '''
    data_type: 'CNV' 'SEG'
    '''
    
    try:
        assert data_type in ['CNV','SEG']
    except:
        print("data type can only be 'CNV' or 'SEG' ")
        exit()
            
    recall_list = []
    precision_list =[]
    cnt = 0 
    for i in p_list:
        eval_dir_list = df_config[df_config['ploidy']==i]['eval_dir'].tolist()
        dc_recall, dc_precision, dc_f1, df_new = collect_data(i,  tool_list, eval_dir_list, data_type)
        if cnt == 0:
            df = df_new
        else:
            df = pd.concat([df, df_new])
        cnt+=1
        recall_list.append(dc_recall)
        precision_list.append(dc_precision)
    df = df.reset_index(drop = True)
    outfile = out_dir + "/"+ data_type + "_stats.csv"
    df.to_csv(outfile, index = False)

    final_info = []
    for p in p_list:
        for tool in tool_list:
            info = stats_one_tool(df, tool,p)
            final_info.extend(info)
    df1 = pd.DataFrame(final_info, columns = ['p','metric','tool','Ncell','mean','std','max','min','range'])
    df1 = df1.sort_values(['p','metric','mean'])
    outfile1 = out_dir + "/"+ data_type + "_stats_agg.csv"
    df1.to_csv(outfile1, index = False)






if not os.path.exists(output_dir):
    os.system("mkdir -p " + output_dir)

df_config = pd.read_csv(input_path)
tool_list = []
p_list = []
for i in range(df_config.shape[0]):
    tool = df_config['tool'][i]
    p= df_config['ploidy'][i]
    if tool not in tool_list:
        tool_list.append(tool)
    if p not in p_list:
        p_list.append(p)


collect('CNV', tool_list, p_list, df_config, output_dir)
collect('SEG', tool_list, p_list, df_config, output_dir)