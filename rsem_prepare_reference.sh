#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=20gb,walltime=24:00:00
#PBS -m abe
#PBS -N RSEM_prepare_reference


module load bowtie/1.0.0
export PATH=$PATH:~/rsem-1.2.7/

cd ${PBS_O_WORKDIR}
#rsem-prepare-reference --transcript-to-gene-map knownIsoforms.txt transcripts.fa.clean.nr.rsem transcripts-rsem
rsem-prepare-reference --transcript-to-gene-map $knownIsoforms $input $output
