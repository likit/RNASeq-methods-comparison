#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=12gb,walltime=12:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Trimming_single_${PBS_JOBID}

cd ${PBS_O_WORKDIR}
perl ~/condetri_v2.1.pl -fastq1=${input} -sc=33 -cutfirst 10
