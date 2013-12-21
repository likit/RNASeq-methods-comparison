#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=64gb,walltime=24:00:00
#PBS -m abe
#PBS -N Velveth_global

cd ${PBS_O_WORKDIR}
velveth global 21,33,2 -fastq -shortPaired ${pe_input} -short ${se_input}
