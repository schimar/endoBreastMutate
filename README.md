# Generic snakemake workflow

======================================================

Snakemake workflow for DNA analysis (on the mach2 HPC cluster with the PBS-torque batch job submission system). This repository was initially based on the [Snakemake Cluster Tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial.git) and the [Software Carpentry lesson repository](https://hpc-carpentry.github.io/hpc-python/17-cluster/).  


======================================================

## conda and other [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s6.yaml)   

create environment from yaml file (in envs/):
```
# run these two once, to create the environment:
conda init bash

module load Anaconda3/2021.04/miniconda-base-2021.04 
source $UIBK_CONDA_PROFILE 

mamba env create -f envs/s6.yaml

# with this, you can activate the environment with all [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml):
conda activate smk6

# (also, when ssh'ing onto mach2, you can activate the env and then do a dry-run of your workflow) 


# if you've added new software to install to the conda environment, then you can update:
conda env update --name ta --file envs/s6.yaml
```
NOTE: the update command didn't quite work in my case, but I had to either ```mamba remove env smk6``` and create the env anew, or use the ```mamba env create``` coimmand again. Do whatever works for you, I suppose...

## Submit the main snakemake job:
```
qsub code/clusterSnakemake.pbs
```




