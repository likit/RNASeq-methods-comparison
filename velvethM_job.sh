#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=256gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_global_merge

cd ${PBS_O_WORKDIR}
velveth global_merged 27 -long global_*/transcripts.fa
