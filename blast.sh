#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N BLAST_${PBS_JOBID}

module load BLAST
cd ${PBS_O_WORKDIR} 

export BLASTDB=/mnt/ls12/preeyanon/blastdb

blastall -p ${program} -b 20 -v 20  -e 0.001 -m 7 -d nr -i ${input} -o ${input}.xml
