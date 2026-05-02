import pandas as pd

tool_list = ['CHISEL', 'HMMcopy', 'rcCAE', 'Ginkgo', 'AneuFinder', 'SeCNV', 'SCOPE', 'CNVeil']


tool = tool_list [ 0]
data = []
for tool in tool_list:
    with open(f"acgh_eval/eval_t10_{tool}.txt",'r') as f:
        s = f.readlines()


    for line in s:
        tp = line.split('_')[0]
        val = [[eval(x),tp,tool] for x in line[:-1].split(',')[1:]]
        data.extend(val)
data = pd.DataFrame(data, columns = ['MSE','subclone','Tool'])

tool_list = ['CHISEL', 'HMMcopy', 'rcCAE', 'Ginkgo', 'AneuFinder', 'SeCNV', 'SCOPE', 'CNVeil']
colors = ['#7f7f7f', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#1f77b4' ,'#e377c2']

color_dict = dict(zip(tool_list, colors))

data['Tool_Color'] = data['Tool'].map(color_dict)

# import seaborn library 
import seaborn as sns 
import matplotlib.pyplot as plt
# load the dataset 
# data = sns.load_dataset('tips') 
# Set the font size for the entire plot
plt.rcParams.update({'font.size': 18})  # Adjust the font size as needed
plt.rcParams['pdf.fonttype'] = 42
# view the dataset 
# print(data.head(5))

# plt.show()

# # g = sns.FacetGrid(data, col='subclone', height=4, aspect=1.5)

# # Map a boxplot to each subplot using sns.boxplot
# # g.map(sns.boxplot, 'tool', 'mse', 'tool',width=0.5, dodge=True)
# tp = 'D'
# sns.set(font_scale=2)
# data1 = data[data['subclone']==tp]
# ax = sns.boxplot(x = data1['subclone'], 
#             y = data1['mse'], 
#             hue = data1['tool'], width = 0.5, dodge=True
#                 )
# ax.set_ylim(0,0.15)
# # ax.get_legend().set_visible(False)

plt.figure(figsize=(15,10))
plt.grid()
ax = sns.boxplot(x = data['subclone'], 
            y = data['MSE'], 
            hue = data['Tool'], width = 0.5, dodge=True,
            palette=color_dict, 
                )

ax.set_ylim(0,4)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f"acgh_bycell.pdf", bbox_inches='tight')
plt.close()


plt.figure(figsize=(15,10))
plt.grid()
ax = sns.boxplot(x = data['subclone'], 
            y = data['MSE'], 
            hue = data['Tool'], width = 0.5, dodge=True,
            palette=color_dict, 
                )
# ax.set_ylim(0,4)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f"acgh_bycell_full.pdf", bbox_inches='tight')


data = data[data['Tool']!='CHISEL']
plt.figure(figsize=(15,10))
plt.grid()
ax = sns.boxplot(x = data['subclone'], 
            y = data['MSE'], 
            hue = data['Tool'], width = 0.5, dodge=True,
            palette=color_dict, 
                )
ax.set_ylim(0,6)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f"acgh_bycell_no_chisel.pdf", bbox_inches='tight')
