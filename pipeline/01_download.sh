#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 16 -n 1 --mem 16gb --out logs/download_sra.%a.log

module load parallel-fastq-dump
module load workspace/scratch

CPU=2
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

FOLDER=sra_data
mkdir -p $FOLDER

N=1
if [ ! -z $SLURM_ARRAY_TASK_ID ]; then
    N=$SLURM_ARRAY_TASK_ID
fi
SRA=$(sed -n ${N}p sralist.tsv);
parallel-fastq-dump -O $FOLDER  --threads $CPU --split-files --gzip --sra-id $SRA
