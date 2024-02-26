# Orginal Vignettes
The package paths are:

```bash
https://github.com/raphael-group/chisel](https://github.com/deepomicslab/SeCNV)
```

# Environments
SeCNV is a Python based method. I installed `conda` and set up the package environment in my path. You can use that environment by adding the code chunk below in your bash before using any functions in SeCNV package:

```bash
ml Anaconda3/5.0.1
source activate /home/liy109/.conda/envs/SeCNV

ml GCC/10.2.0
ml SAMtools/1.12
```

# For scDNA-seq data
### Prepare the data for SeCNV
To run SeCNV, the bigwig files, hg19_mappability.bigWig and hg38_mappability.bigWig, should be downloaded from (https://drive.google.com/drive/folders/1XGuXa9muRtOiAjtz1J4gWO5Qk3c5-PwS) and put under the Script folder. And picard.jar file should also under the dataset folder.

You could see the code chunk for generation of KTN126 in the preparation of the data.

### Use SeCNV to infer CNV
After preparation, we can use tool `SeCNV` to infer CNV:

```bash
module load Anaconda3/5.0.1
source /home/liy109/.conda/envs/SeCNV

ml GCC/10.2.0
ml SAMtools/1.12

cd Scripts 
python SeCNV.py /data/maiziezhou_lab/Datasets/singlecell_data/SeCNV/T10_hg19/Input /data/maiziezhou_lab/Datasets/singlecell_data/SeCNV/T10_hg19/output hg19.fa
```

Note that, please make sure the genome reference is the same during the whole procedure. For different datasets, we only have to change the path of input folder, output folder and reference file.
