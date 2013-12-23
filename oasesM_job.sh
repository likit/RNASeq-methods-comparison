#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=256gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N OasesM_${PBS_JOBID}

cd ${PBS_O_WORKDIR}
~/oases_0.2.06/oases global_merged -merge yes
