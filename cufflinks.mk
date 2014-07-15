run-quality-trim-pe:

	qsub -v left="reads/line7u.pe.1,right=reads/line7u.pe.2" \
		$(protocol)/quality_trim_pe_job.sh
	qsub -v left="reads/line7i.pe.1,right=reads/line7i.pe.2" \
		$(protocol)/quality_trim_pe_job.sh
	qsub -v left="reads/line6u.pe.1,right=reads/line6u.pe.2" \
		$(protocol)/quality_trim_pe_job.sh
	qsub -v left="reads/line6i.pe.1,right=reads/line6i.pe.2" \
		$(protocol)/quality_trim_pe_job.sh

run-quality-trim-se:

	for r in reads/*.se.fq; do qsub -v input="$$r" \
		$(protocol)/quality_trim_se_job.sh; done

run-tophat-se:

	# requires Tophat 2.0.9 and Bowtie 2.1.0
	cd tophat; qsub -v "outdir=line6u_se,index=gal4selected,\
		input=../reads/line6u.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v "outdir=line6i_se,index=gal4selected,\
		input=../reads/line6i.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v "outdir=line7u_se,index=gal4selected,\
		input=../reads/line7u.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v "outdir=line7i_se,index=gal4selected,\
		input=../reads/line7i.fq_trim.fastq" $(protocol)/tophat_se_job.sh

run-tophat-pe:

	# requires Tophat 2.0.9 and Bowtie 2.1.0
	cd tophat; qsub -v "outdir=line6u_pe,index=gal4selected,\
		left=../reads/line6u.1_trim1.fastq,right=../reads/line6u.1_trim2.fastq,\
		unpaired=../reads/line6u.1_trim_unpaired.fastq" \
		$(protocol)/tophat_pe_job.sh 

	cd tophat; qsub -v "outdir=line6i_pe,index=gal4selected,\
		left=../reads/line6i.1_trim1.fastq,right=../reads/line6i.1_trim2.fastq,\
		unpaired=../reads/line6i.1_trim_unpaired.fastq" \
		$(protocol)/tophat_pe_job.sh 

	cd tophat; qsub -v "outdir=line7u_pe,index=gal4selected,\
		left=../reads/line7u.1_trim1.fastq,right=../reads/line7u.1_trim2.fastq,\
		unpaired=../reads/line7u.1_trim_unpaired.fastq" \
		$(protocol)/tophat_pe_job.sh 

	cd tophat; qsub -v "outdir=line7i_pe,index=gal4selected,\
		left=../reads/line7i.1_trim1.fastq,right=../reads/line7i.1_trim2.fastq,\
		unpaired=../reads/line7i.1_trim_unpaired.fastq" \
		$(protocol)/tophat_pe_job.sh 

run-cufflinks:

	# requires Cufflinks 2.1.1
	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		$(protocol)/cufflinks_job.sh; echo $$d; done

run-cuffmerge-ref:

	# requires Cufflinks 2.1.1
	#../Gallus_UCSC_ensembl_73.gtf.removed
	cd tophat; cuffmerge -g ../Gallus_gallus.Galgal4.73.removed.gtf
		-o merged_cuff_ref -s gal4selected.fa -p 4 $(protocol)/merge_list.txt
