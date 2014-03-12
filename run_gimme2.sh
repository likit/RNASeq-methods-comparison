#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=32gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Gimme_${PBS_JOBID}

module load bxPython
module load pygr
module load matplotlib

cd ${PBS_O_WORKDIR}
python ~/gimme/src/gimme.py -r ${ref} ${input1} ${input2} > ${input1}.merged.bed 2>${input}.merged.gimme.log
