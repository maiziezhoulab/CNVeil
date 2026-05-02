import os
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define the parent directory where tool subfolders are located
parent_folder = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/T10_EVALUATION/Evaluation"
# out_dir = "./eval_t10/"
out_dir = "./evaluation/eval_T10/"

if not os.path.exists(out_dir):
    os.system("mkdir -p " + out_dir)
# Initialize lists to store data
data = []

# Loop through each tool's folder
for tool in os.listdir(parent_folder):
    if tool!='CNVeil':
        tool_path = os.path.join(parent_folder, tool)
    else:
        tool_path = "./T10/eval_acgh"
    print(tool_path)


    # Ensure it's a directory
    if os.path.isdir(tool_path):
        mse_file = os.path.join(tool_path, "MSE.pkl")
        pearson_file = os.path.join(tool_path, "pearson_coef.pkl")

        # Load MSE data
        if os.path.exists(mse_file):
            with open(mse_file, "rb") as f:
                mse_dict = pickle.load(f)
                for subclone, values in mse_dict.items():
                    for value in values:
                        data.append({"Subclone": subclone, "Metric": "MSE", "Value": value, "Tool": tool})

        # Load Pearson data
        if os.path.exists(pearson_file):
            with open(pearson_file, "rb") as f:
                pearson_dict = pickle.load(f)
                for subclone, values in pearson_dict.items():
                    for value in values:
                        data.append({"Subclone": subclone, "Metric": "Pearson", "Value": value, "Tool": tool})

# Convert list to DataFrame
df = pd.DataFrame(data)


# print(df)
# exit()
# Create the subplots for MSE and Pearson
tool_list = [ 'HMMcopy', 'Chisel', 'Ginkgo', 'AneuFinder', 'SeCNV', 'SCOPE', 'CNVeil']
colors = [ '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#1f77b4' ,'#e377c2']

color_dict = dict(zip(tool_list, colors))

fig, axes = plt.subplots(2, 1, figsize=(14, 14))

# Boxplot for MSE
sns.boxplot(x="Subclone", y="Value", hue="Tool", data=df[df["Metric"] == "MSE"], ax=axes[0], palette=color_dict)
axes[0].set_title("MSE Comparison")
axes[0].set_ylabel("MSE Value")
axes[0].set_xlabel("Subclone")

# Boxplot for Pearson
sns.boxplot(x="Subclone", y="Value", hue="Tool", data=df[df["Metric"] == "Pearson"], ax=axes[1], palette=color_dict)
axes[1].set_title("Pearson Coefficient Comparison")
axes[1].set_ylabel("Pearson Coefficient Value")
axes[1].set_xlabel("Subclone")

# Adjust legend position
axes[0].legend(title="Tool", bbox_to_anchor=(1.05, 1), loc='upper left')
axes[1].legend(title="Tool", bbox_to_anchor=(1.05, 1), loc='upper left')

# Show the plot
plt.tight_layout()
plt.savefig(f"{out_dir}/eval_T10_acgh_boxplot.pdf")


df1 = df.groupby(['Metric', 'Tool','Subclone'])['Value'].mean().reset_index()
df1 = df1.sort_values(['Metric','Subclone','Value']).reset_index(drop= True)[['Metric','Subclone','Value','Tool']]

df2 = df.groupby(['Metric', 'Tool',])['Value'].mean().reset_index()
df2 = df2.sort_values(['Metric','Value']).reset_index(drop= True)[['Metric','Value','Tool']]

# print(df1)
df1.to_csv(f"{out_dir}/eval_T10_acgh_ranking_subclone.tsv", sep = '\t')
df2.to_csv(f"{out_dir}/eval_T10_acgh_ranking_overall.tsv", sep = '\t')

# Write to Excel with two sheets
with pd.ExcelWriter(f"{out_dir}/eval_T10_acgh_ranking.xlsx") as writer:
    df1.to_excel(writer, sheet_name='Subclone', index=False)
    df2.to_excel(writer, sheet_name='Overall', index=False)
# df1.to_excel(f"{out_dir}/eval_T10_acgh_ranking.xlsx", index=False)  # Save without index
