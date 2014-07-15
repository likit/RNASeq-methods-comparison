interleave-reads:

	cd assembly; cat ../reads/*trim1.fastq >> pe1.fastq
	cd assembly; cat ../reads/*trim2.fastq >> pe2.fastq
	cd assembly; cat ../reads/*trim_unpaired.fastq \
		../reads/*fq_trim.fastq >> single.fastq
	cd assembly; $(protocol)/shuffleSequences_fastq.pl \
		pe1.fastq pe2.fastq paired.fastq

run-velveth:

	# requires Velvet 1.2.03
	cd assembly; qsub -v \
		pe_input="paired.fastq",se_input="single.fastq" \
		$(protocol)/velveth_job.sh

run-velvetg:

	# requires Velvet 1.2.03
	cd assembly; qsub $(protocol)/velvetg_job.sh

run-oases:

	# requires Oases 0.2.06
	cd assembly; qsub $(protocol)/oases_job.sh
