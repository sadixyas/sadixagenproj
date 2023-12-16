#!/usr/bin/bash -l
##to download mrna.fasta of reference genome of our species of interest
curl -o cds.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/507/305/GCA_000507305.1_ASM50730v1/GCA_000507305.1_ASM50730v1_genomic.fna.gz

zmore cds.fna.gz

##to download input SRA using loop

nano 00_download.sh
module load sratoolkit
##config https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration


