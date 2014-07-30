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

	cd tophat; \
		cuffmerge -g ../Gallus_gallus.Galgal4.73.removed.gtf \
		-o merged_cuff_ref -s ../gal4selected.fa -p 4 $(protocol)/merge_list.txt

annotate-cufflinks:

	cd tophat/merged_cuff_ref; \
		python $(protocol)/gene-rep-cufflinks.py \
		merged-ref.transcripts.fa > merged-ref-genes.fa

	cd tophat/merged_cuff_ref; \
		qsub -v "db=Gallus_prot,input=merged-ref-genes.fa,\
		program=blastx,output=merged-ref-genes-gga.xml" \
		$(protocol)/blast.sh

	cd tophat/merged_cuff_ref; \
		qsub -v "db=Human_prot,input=merged-ref-genes.fa,\
		program=blastx,output=merged-ref-genes-hsa.xml" \
		$(protocol)/blast.sh

prepare-reference-cufflinks:

	# only genes with corresponding Ensembl
	cd tophat/merged_cuff_ref; \
		python $(protocol)/get_top_hits.py merged-ref-genes-gga.xml > merged-ref-gga-tophits.txt

	cd tophat/merged_cuff_ref; \
		cat merged.rsem.gtf | python $(protocol)/cufflinks-to-known-isoforms.py > knownIsoforms.txt

	cd tophat/merged_cuff_ref; \
		python $(protocol)/get_best_ensembl_hits_cufflinks.py \
		merged-ref-gga-tophits.txt merged-ref-genes.fa \
		knownIsoforms.txt > cufflinks-gga.fa

	cd tophat/merged_cuff_ref; \
		qsub -v "input=cufflinks-gga.fa,knownIsoforms=knownIsoforms.txt,\
		output=cufflinks-gga-rsem" $(protocol)/rsem_prepare_reference.sh

run-rsem-calc-expression-cufflinks:

	# only genes with corresponding Ensembl
	cd tophat/merged_cuff_ref; \
		qsub -v "index=cufflinks-gga-rsem,input_read=../../reads/line7u.se.fq,\
		sample_name=line7u-single-rsem" $(protocol)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \
		qsub -v "index=cufflinks-gga-rsem,input_read=../../reads/line7i.se.fq,\
		sample_name=line7i-single-rsem" $(protocol)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \
		qsub -v "index=cufflinks-gga-rsem,input_read1=../../reads/line7u.pe.1,\
		input_read2=../../reads/line7u.pe.2,sample_name=line7u-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

	cd tophat/merged_cuff_ref; \
		qsub -v "index=cufflinks-gga-rsem,input_read1=../../reads/line7i.pe.1,\
		input_read2=../../reads/line7i.pe.2,sample_name=line7i-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-cufflinks:

	cd tophat/merged_cuff_ref; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results \
		line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix

	cd tophat/merged_cuff_ref; rsem-run-ebseq \
		line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes

	cd tophat/merged_cuff_ref; rsem-control-fdr \
		line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-tophits-degenes:

	cd tophat/merged_cuff_ref; \
	python $(protocol)/tophits-to-degenes-cufflinks.py \
		line7u_vs_i.degenes.fdr.05 merged-ref-gga-tophits.txt \
		knownIsoforms.txt > line7u_vs_i.degenes.fdr.05.gga.tophits

	cd tophat/merged_cuff_ref; \
		python $(protocol)/get_top_hits.py merged-ref-genes-hsa.xml > merged-ref-hsa-tophits.txt

	cd tophat/merged_cuff_ref; \
	python $(protocol)/tophits-to-degenes-cufflinks.py \
		line7u_vs_i.degenes.fdr.05 merged-ref-hsa-tophits.txt \
		knownIsoforms.txt > line7u_vs_i.degenes.fdr.05.hsa.tophits

