#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Velveth_local_${PBS_JOBID}
#PBS -A ged-intel11

cd ${PBS_O_WORKDIR}
~/velvet_1.2.03/velveth ${outdir} 21,33,2 -bam -shortPaired ${input}
