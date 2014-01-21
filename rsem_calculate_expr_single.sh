#!/bin/sh -login
#PBS -l nodes=1:ppn=8,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N RSEM_calc_expr_${sample_name}

cd ${PBS_O_WORKDIR}
module load bowtie
/mnt/home/preeyano/rsem-1.2.7/rsem-calculate-expression --time -p 8 ${input_read} ${index} ${sample_name}
