#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 16  --out logs/kallisto.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

module load kallisto
INDIR=fastq
OUTDIR=results/kallisto
REFERENCE_GENOME="/bigdata/stajichlab/sadikshs/gen220/genproject/reference_genome/ncbi_dataset/data/GCA_000507305.1/GCA_000507305.1_ASM50730v1_genomic.fna"
FASTQ_FILES=("SRR12010079.fastq.gz" "SRR12010080.fastq.gz" "SRR12010081.fastq.gz" "SRR12010082.fastq.gz")
SPECIES=Algae
INPUT=input
source config.txt

IDX=db/${SPECIES}.idx
TX=db/${SPECIES}_mRNA.fasta
SAMPLEFILE=samples.csv

mkdir -p $OUTDIR
if [ ! -f $IDX ]; then
    kallisto index -i $IDX $REFERENCE_GENOME
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE CONDITION REPLICATION READBASE
do
 OUTNAME=$CONDITION.r${REP}
 if [ ! -f $OUTDIR/$OUTNAME/abundance.h5 ]; then
     kallisto quant -i $IDX -o $OUTDIR/$OUTNAME -t $CPU --bias $INDIR/${READBASE}_1.fastq.gz $INDIR/${READBASE}_2.fastq.gz
 fi
done

