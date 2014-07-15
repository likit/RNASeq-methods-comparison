interleave-reads:

	cd assembly; cat ../reads/*trim1.fastq >> pe1.fastq
	cd assembly; cat ../reads/*trim2.fastq >> pe2.fastq
	cd assembly; cat ../reads/*trim_unpaired.fastq \
		../reads/*fq_trim.fastq >> single.fastq
	cd assembly; $(protocol)/shuffleSequences_fastq.pl \
		pe1.fastq pe2.fastq paired.fastq

run-velveth:

	cd assembly; qsub -v \
		pe_input="paired.fastq",se_input="single.fastq" \
		$(protocol)/velveth_job.sh

run-velvetg:

	cd assembly; qsub $(protocol)/velvetg_job.sh

run-oases:

	cd assembly; qsub $(protocol)/oases_job.sh
