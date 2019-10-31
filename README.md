# CutAndRun
Snakemake process for trimming and alignment of Cut&amp;Run sequencing

# For use in linux
# Developed with snakemake (v3.13.3)

##Create the conda environment 
conda env create -f environment_cutrun.yaml

##Install ngsutilsj
https://github.com/compgen-io/ngsutilsj
##Note: rule sort_combined_lengths will require a correct path for ngsutilsj

##Create directory structure
#Within the working directory:
mkdir ./fastq/
#Place sample fastq.gz files in the ./fastq directory

