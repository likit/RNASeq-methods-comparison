#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N BLAST2GO_${PBS_JOBID}

module load blast2go

cd ${PBS_O_WORKDIR} 

java es.blast2go.prog.B2GAnnotPipe -goslim -annex -in ${input} -ips interpro/ -dat -annot -prop b2gPipe.properties -out ${outdir} -img
