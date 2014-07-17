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

run-oasesM-velveth:

	cd assembly; qsub $(protocol)/velvethM_job.sh

run-oasesM-velvetg:

	cd assembly; qsub $(protocol)/velvetgM_job.sh

run-oasesM:

	cd assembly; qsub $(protocol)/oasesM_job.sh

clean-transcripts:

	# -A needed to keep poly-A tail
	cd assembly/global_merged; ~/seqclean-x86_64/seqclean transcripts.fa -c 8 -A -o transcripts.fa.clean
	qsub -v input="assembly/global_merged/transcripts.fa.clean",output="assembly/global_merged/transcripts.fa.clean.nr",c="1.0" \
		$(protocol)/cdhit_job.sh

remove-redundant-seq:

	cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh
