#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=48gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Velvetg_local_${PBS_JOBID}
#PBS -A ged-intel11

cd ${PBS_O_WORKDIR}
~/velvet_1.2.03/velvetg ${indir} -clean yes -read_trkg yes -ins_length 175
