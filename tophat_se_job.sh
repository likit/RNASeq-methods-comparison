#!/bin/sh -login
#PBS -l nodes=1:ppn=4,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Tophat_single_${PBS_JOBID}

module load bowtie2/2.1.0
cd ${PBS_O_WORKDIR}

~/tophat-2.0.9.Linux_x86_64/tophat -p 4 -o ${outdir} ${index} ${input}
