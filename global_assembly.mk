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
	cd assembly/global_merged; seqclean transcripts.fa -c 8 -A -o transcripts.fa.clean
	qsub -v "input=assembly/global_merged/transcripts.fa.clean,\
		output=assembly/global_merged/transcripts.fa.clean.nr,c=1.0" \
		$(protocol)/cdhit_job.sh

annotate-assembly:

	cd assembly/global_merged; python $(protocol)/gene-rep-velvet.py \
		transcripts.fa.clean.nr > assembly-genes.fa
	cd assembly/global_merged; \
		qsub -v "db=Gallus_prot,input=assembly-genes.fa,program=blastx,\
		output=assembly-genes-gga.xml" $(protocol)/blast.sh
	cd assembly/global_merged; \
		qsub -v "db=Human_prot,input=assembly-genes.fa,program=blastx,\
		output=assembly-genes-hsa.xml" $(protocol)/blast.sh

rsem-prepare-reference-assembly:

	cd assembly/global_merged; \
		python $(protocol)/get_top_hits.py assembly-genes-gga.xml > assembly-gga-tophits.txt

	cd assembly/global_merged; \
		python $(protocol)/get_best_ensembl_hits_assembly.py assembly-gga-tophits.txt \
		transcripts.fa.clean.nr > transcripts-gga.fa

	cd assembly/global_merged; \
		cat transcripts-gga.fa | python $(protocol)/prepare-transcripts.py \
		transcripts-gga-rsem.fa knownIsoforms.gga.txt

	cd assembly/global_merged; \
		qsub -v "input=transcripts-gga-rsem.fa,knownIsoforms=knownIsoforms.gga.txt,\
		output=transcripts-gga-rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-expression-assembly:

	cd assembly/global_merged; \
		qsub -v "input_read=../../reads/line7u.se.fq,\
		sample_name=line7u-single-rsem,\
		index=transcripts-gga-rsem" $(protocol)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; \
		qsub -v "input_read=../../reads/line7i.se.fq,\
		sample_name=line7i-single-rsem,\
		index=transcripts-gga-rsem" $(protocol)/rsem_calculate_expr_single.sh

	cd assembly/global_merged; \
		qsub -v "input_read1=../../reads/line7u.pe.1,\
		input_read2=../../reads/line7u.pe.2,\
		sample_name=line7u-paired-rsem,\
		index=transcripts-gga-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; \
		qsub -v "input_read1=../../reads/line7i.pe.1,\
		input_read2=../../reads/line7i.pe.2,\
		sample_name=line7i-paired-rsem,\
		index=transcripts-gga-rsem" $(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-assembly:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results  \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results  \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix

	cd assembly/global_merged; \
		rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 \
		line7u_vs_i.degenes

	cd assembly/global_merged; \
		rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-tophits-degenes:

	cd assembly/global_merged; \
		python $(protocol)/tophits-to-degenes-assembly.py \
		line7u_vs_i.degenes.fdr.05 assembly-gga-tophits.txt \
		> line7u_vs_i.degenes.fdr.05.gga.tophits

	cd assembly/global_merged; \
		python $(protocol)/get_top_hits.py assembly-genes-hsa.xml \
		> assembly-hsa-tophits.txt

	cd assembly/global_merged; \
		python $(protocol)/tophits-to-degenes-assembly.py \
		line7u_vs_i.degenes.fdr.05 assembly-hsa-tophits.txt \
		> line7u_vs_i.degenes.fdr.05.hsa.tophits

