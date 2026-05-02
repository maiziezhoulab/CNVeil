import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import os
def read_matrix(file_name):
    matrix = []
    chr_name = []
    bin_list = []
    count = 0

    with open(file_name, 'r') as f:
        for line in f:
            count += 1
            line = line.strip().split('\t')

            if count > 1:
                matrix.append(list(map(int, line[3:])))
                chr_name.append(line[0])
                bin_list.append(':'.join(line[:3]))
            else:
                sample_list = line[3:]

    return matrix, chr_name, bin_list, sample_list

def get_gini(Y):
    n_cols = Y.shape[1]
    Gini = np.empty(n_cols)

    for i in range(n_cols):
        y = np.sort(Y[:, i])
        x = np.concatenate(([0], np.arange(1, len(y) + 1) / len(y)))
        y = np.concatenate(([0], np.cumsum(y) / np.sum(y)))
        Gini[i] = 2 * round(0.5 - np.trapz(y, x), 4)

    return Gini


def get_norm_cell(data_set, normal_cell_file=None, quan_vec=0.4, gini_threshold=0.12):
    # data_set = np.transpose(data_set)
    # print(data_set.shape)
    data_set = np.array(data_set)
    
    pca = PCA(n_components=1, whiten=True)
    
    pca_result = pca.fit_transform(np.array(data_set).T)
    # print(pca.components_.shape)
    


    var = pca.explained_variance_
    # print(var)
    rotation = pca.components_.flatten()
    # print(rotation)
    # print(sum(rotation))

    sorted_bins = np.argsort(np.abs(rotation))[::-1]
    # print(sorted_bins)
    bin_index = sorted_bins[ rotation[sorted_bins] >= np.quantile(rotation, quan_vec)]
    # print(len(bin_index))
    # print(sorted(bin_index)[:100])
    # print(data_set.shape)
    # bin_index = [1,2]
    Y_selected = data_set[ bin_index,:]
    

    gini_coefficients = get_gini(Y_selected)
    # print(gini_coefficients)

    if normal_cell_file is None:
        norm_index = np.where(gini_coefficients < gini_threshold)[0]
        abnorm_index = np.where(gini_coefficients >= gini_threshold)[0]
        # print(sorted(norm_index+1))
        print(f"Total cells: {data_set.shape[1]} \nNumber of normal cells: {len(norm_index)}")
    else:
        print(f"Using predefined normal cells from {normal_cell_file}")
        predefined_normal_cells = pd.read_table(normal_cell_file, header=None, names=["sample"], delimiter="\t")
        predefined_normal_cells = predefined_normal_cells["sample"].tolist()
        norm_index = np.where(np.isin(data_set.columns, predefined_normal_cells))[0]
        abnorm_index = np.where(~np.isin(data_set.columns, predefined_normal_cells))[0]

    return {"norm_index": norm_index, "abnorm_index": abnorm_index}



def Bias_norm(Y, norm_cell_index=None):
    Y = np.transpose(Y)

    if norm_cell_index is not None and len(norm_cell_index) > 0:
        norm_cell_Y = Y[norm_cell_index, :]
        bias_matrix = np.zeros((norm_cell_Y.shape[0], norm_cell_Y.shape[1]))

        for i in range(norm_cell_Y.shape[0]):
            cell = norm_cell_Y[i, :]
            median = np.median(cell)
            bias_list = cell / median
            bias_matrix[i, :] = bias_list

        ave_bias = np.mean(bias_matrix, axis=0)
        ave_bias[ave_bias == 0] = 1
    else:
        with open("bias.txt", "r") as file:
            ave_bias = np.array([float(value) for value in file.read().split()])

    gc_nor_Y = Y / np.tile(ave_bias, (Y.shape[0], 1))
    return np.transpose(gc_nor_Y)

def normal_cell_correction(work_path):
    cov_bed_path = os.path.join(work_path, "genome_cov.bed")
    normal_cell_file = None
    quan_vec = 0.4
    gini_threshold = 0.12

    # Read data and perform the operations
    read_data = read_matrix(cov_bed_path)
    Y = read_data[0]
    chr_name = read_data[1]
    bin_list = read_data[2]
    sample_list = read_data[3]

    var_norm_cell = get_norm_cell(Y, normal_cell_file, quan_vec, gini_threshold)
    cov_matrix = Bias_norm(Y, var_norm_cell["norm_index"])

    # Create a DataFrame from the cov_matrix
    df_data = pd.DataFrame(cov_matrix)

    # Save the DataFrame as a CSV file
    df_data.to_csv(os.path.join(work_path,"RDmatrix.csv"), index=False)


if __name__ == '__main__':
    # Define the paths and parameters
    work_path = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/pipeline_test/Real_Data/T16/"
    normal_cell_correction(work_path)









