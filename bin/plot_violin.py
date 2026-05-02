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
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns





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

def collect_data(p, tool_list, eval_dir_list, data_type):
    # tool_list = [
    #              "ginkgo",
    #              "hmmcopy",
    #              "scope",
    #              "secnv"]
    if data_type == 'CNV':
        try:
            file_list = [ eval_dir_list[i]+'/stats_by_cell_2000000_matchcn.csv' for i in range(len(tool_list))]
        except:
            print(len(eval_dir_list))
            print(len(tool_list))
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
#     return df




def plot_violin(data_type, tool_list, p_list, df_config, out_dir):
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
        tool_list_p = df_config[df_config['ploidy']==i]['tool'].tolist()
        print(tool_list_p)
        dc_recall, dc_precision, dc_f1, df_new = collect_data(i,  tool_list_p, eval_dir_list, data_type)
        if cnt == 0:
            df = df_new
        else:
            df = pd.concat([df, df_new])
        cnt+=1
        recall_list.append(dc_recall)
        precision_list.append(dc_precision)
    df = df.reset_index(drop = True)


    # with the dataframe, now draw the violin 
    rc={'font.size': 24, 'axes.labelsize': 16, 'legend.fontsize': 24.0,
    'axes.titlesize': 28, 'xtick.labelsize': 16, 'ytick.labelsize': 16}
    plt.rcParams['pdf.fonttype'] = 42
    sns.set(rc=rc)
    sns.set(context="paper", style="ticks")
    g = sns.FacetGrid(df, row="metric", sharey=False,  aspect=4)
    # print(df)
    g = g.map(sns.violinplot, "Ploidy", "value", "Tool", cut=0, inner=None, scale="width", split=False, palette=color_dict , saturation=1).despine(left=True)
    # Set axis labels & ticks #
    fontsize=14
    fontsize1=14
    fontsize2=14
    g.fig.get_axes()[0].set_xticklabels(p_list, size=fontsize)
    g.fig.get_axes()[0].set_xlabel("")
    g.fig.get_axes()[0].set_yticks(np.round(np.arange(0,1,0.1),1))
    g.fig.get_axes()[0].set_yticklabels(np.round(np.arange(0,1,0.1),1), size=fontsize)
    g.fig.get_axes()[0].set_ylabel("Recall", size=fontsize2)
    g.fig.get_axes()[0].spines["left"].set_visible(True)
    
    g.fig.get_axes()[1].set_xticklabels(p_list, size=fontsize)
    g.fig.get_axes()[1].set_xlabel("")
    g.fig.get_axes()[1].set_yticks(np.round(np.arange(0,1,0.1),1))
    g.fig.get_axes()[1].set_yticklabels(np.round(np.arange(0,1,0.1),1), size=fontsize)
    g.fig.get_axes()[1].set_ylabel("Precision", size=fontsize2)
    g.fig.get_axes()[1].spines["left"].set_visible(True)

    g.fig.get_axes()[2].set_xticklabels(p_list, size=fontsize)
    g.fig.get_axes()[2].set_xlabel("")
    g.fig.get_axes()[2].set_yticks(np.round(np.arange(0,1,0.1),1))
    g.fig.get_axes()[2].set_yticklabels(np.round(np.arange(0,1,0.1),1), size=fontsize)
    g.fig.get_axes()[2].set_ylabel("F1", size=fontsize2)
    g.fig.get_axes()[2].spines["left"].set_visible(True)


    # Set legend #

    handles, labels = g.fig.get_axes()[0].get_legend_handles_labels()

    g.fig.get_axes()[0].legend(handles, labels, loc='upper left', fontsize='x-large', bbox_to_anchor = (1,1))
    # Fixing titles #

    #----------------- uncomment this if you do not want subfigure title
    g.fig.get_axes()[0].set_title("")
    g.fig.get_axes()[1].set_title("")
    g.fig.get_axes()[2].set_title("")
    g.set_titles(size = 16)
    # plt.title(f"{data_type} evaluation", size = 18)
    plt.xlabel("Simulation ploidy", size = fontsize2 )

    fig_file = out_dir + f"/combined_quan_{data_type}.png"
    fig_file1 = out_dir + f"/combined_quan_{data_type}.pdf"
    
    plt.tight_layout()
    plt.savefig(fig_file, dpi=400)
    plt.savefig(fig_file1, dpi=400)
    plt.show()
    plt.close()
    
    return 

tool_list = [ 'HMMcopy', 'rcCAE', 'Ginkgo', 'AneuFinder', 'SeCNV', 'SCOPE', 'CNVeil']
colors = [ '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#1f77b4' ,'#e377c2']

color_dict = dict(zip(tool_list, colors))

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

print(tool_list)
plot_violin('CNV', tool_list, p_list, df_config, output_dir)
plot_violin( 'SEG', tool_list, p_list, df_config, output_dir)







