#!/bin/sh -login
#PBS -l nodes=1:ppn=4,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N RSEM_calc_expr_${sample_name}

cd ${PBS_O_WORKDIR}
module load bowtie
/mnt/home/preeyano/rsem-1.2.7/rsem-calculate-expression --paired-end --output-genome-bam --time -p 4 ${input_read1} ${input_read2} galGal4-removed ${sample_name}