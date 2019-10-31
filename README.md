# CutAndRun
Snakemake process for trimming and alignment of Cut&amp;Run sequencing

For use in linux

Developed with snakemake (v3.13.3)

# Create the conda environment 
conda env create -f environment_cutrun.yaml

# Install auxillary programs and create directory structure

Install ngsutilsj
(https://github.com/compgen-io/ngsutilsj)

Note: rule sort_combined_lengths will require a correct path for ngsutilsj

Create directory structure within the working directory:

mkdir ./fastq/

Place sample fastq.gz files in the ./fastq directory
```
.
├── config.yaml
├── cluster.json
├── fastq
│   ├── Sample1_250k_R1_001.fastq.gz
│   ├── Sample1_250k_R1_001.fastq.gz
│   ├── Sample2_250k_R1_001.fastq.gz
│   └── Sample2_250k_R2_001.fastq.gz
├── run_snakemake_cluster.sh
└── Snakefile
```
##Run 

conda activate cutrun ##activate environment

snakemake -nr --snakefile Snakefile ##shows what commands will be run

nohup sh run_snakemake_cluster.sh & ## to run commands via bsub cluster system
