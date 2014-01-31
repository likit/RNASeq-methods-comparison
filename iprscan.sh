#!/bin/sh -login
#PBS -l nodes=1:ppn=5,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N InterproScan_${PBS_JOBID}

module load iprscan
cd ${PBS_O_WORKDIR} 

iprscan5 -i ${input}
