##############################################
############# Ensembl annotations #############
###############################################

rsem-prepare-ensembl-reference:

	python $(protocol)/gtf_to_known_isoforms.py \
		Gallus_gallus.Galgal4.73.removed.gtf > Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt

	rsem-prepare-reference --gtf Gallus_gallus.Galgal4.73.removed.gtf \
		--transcript-to-gene-map Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt \
		gal4selected.fa ensembl_genes

rsem-calc-expression:

	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem",index="ensembl_genes" \
		~/rnaseq-comp-protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem",index="ensembl_genes" \
		~/rnaseq-comp-protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="ensembl_genes" \
		~/rnaseq-comp-protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="ensembl_genes" \
		~/rnaseq-comp-protocol/rsem_calculate_expr_paired.sh

run-ebseq-line7:

	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-ensembl:

	python ~/rnaseq-comp-protocol/rsem-output-ensbl-to-fasta.py line7u_vs_i.degenes.fdr.05 \
		ensembl_genes.transcripts.fa Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt \
		> line7u_vs_i.degenes.fdr.05.fa

run-blast-ensembl-gallus:

	python $(protocol)/gene-rep-ensbl.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)l/blast.sh

run-blast-ensembl-human:

	mkdir Human_blast
	qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)/blast.sh

run-goseq-ensembl-gallus:

	Rscript $(protocol)/goseq_ensembl_gallus.R

###############################################
##### De novo assembly with Velvet+OasesM #####
###############################################

run-quality-trim-pe:

	qsub -v left="reads/line7u.pe.1,right=reads/line7u.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line7i.pe.1,right=reads/line7i.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line6u.pe.1,right=reads/line6u.pe.2" protocol/quality_trim_pe_job.sh
	qsub -v left="reads/line6i.pe.1,right=reads/line6i.pe.2" protocol/quality_trim_pe_job.sh

run-quality-trim-se:

	for r in reads/*.se.fq; do qsub -v input="$$r" protocol/quality_trim_se_job.sh; done

interleave-reads:

	cd assembly; ~/velvet_1.2.03/shuffleSequences_fastq.pl pe.1.fastq pe.2.fastq paired.fastq

run-velveth:

	cd assembly; qsub -v pe_input="paired.fastq",se_input="single.fastq" $(protocol)/velveth_job.sh

run-velvetg:

	cd assembly; qsub $(protocol)/velvetg_job.sh

run-oases:

	cd assembly; qsub $(protocol)/oases_job.sh

run-oasesM:

	cd assembly; qsub $(protocol)/velvethM_job.sh
	cd assembly; qsub $(protocol)/velvetgM_job.sh
	cd assembly; qsub $(protocol)/oasesM_job.sh

clean-transcripts:

	# -A needed to keep poly-A tail
	cd assembly/global_merged; ~/seqclean-x86_64/seqclean transcripts.fa -c 8 -A -o transcripts.fa.clean
	qsub -v input="assembly/global_merged/transcripts.fa.clean",output="assembly/global_merged/transcripts.fa.clean.nr",c="1.0" \
		$(protocol)/cdhit_job.sh

annotate-global-asm:

	# cd assembly/global_merged; python ~/rnaseq-comp-protocol/gene-rep-velvet.py transcripts.fa.clean.nr > genes.fa
	cd assembly/global_merged; \
		qsub -v db="Gallus_prot",input="genes.fa",program="blastx",output="genes.xml" ~/rnaseq-comp-protocol/blast.sh

rsem-prepare-reference-global-asm-ensembl-matched:

	cd assembly/global_merged; \
	python $(protocol)/get_best_ensembl_hits_assembly.py genes.tophits.txt \
	transcripts.fa.clean.nr > transcripts.ensembl-matched.fa

	cd assembly/global_merged; cat transcripts.ensembl-matched.fa | python $(protocol)/prepare-transcripts.py \
		transcripts.ensembl-matched.rsem.fa knownIsoforms.ensembl-matched.txt

	cd assembly/global_merged; qsub -v "input=transcripts.ensembl-matched.rsem.fa,\
		knownIsoforms=knownIsoforms.ensembl-matched.txt,\
		output=transcripts.ensembl-matched-rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-expression-global-asm-ensembl-matched:

	cd assembly/global_merged; qsub -v input_read="../../reads/line7u.se.fq",sample_name="line7u-single-ensembl-matched-rsem",index="transcripts.ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7i.se.fq",sample_name="line7i-single-ensembl-matched-rsem",index="transcripts.ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_single.sh

	cd assembly/global_merged; qsub -v input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-ensembl-matched-rsem",index="transcripts.ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-ensembl-matched-rsem",index="transcripts.ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-global-asm-ensembl-matched:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line7u-single-ensembl-matched-rsem.genes.results  \
		line7u-paired-ensembl-matched-rsem.genes.results line7i-single-ensembl-matched-rsem.genes.results  \
		line7i-paired-ensembl-matched-rsem.genes.results > line7u_vs_i.gene-ensembl-matched.counts.matrix
	cd assembly/global_merged; rsem-run-ebseq line7u_vs_i.gene-ensembl-matched.counts.matrix 2,2 \
		line7u_vs_i.ensembl-matched.degenes
	cd assembly/global_merged; rsem-control-fdr line7u_vs_i.ensembl-matched.degenes 0.05 line7u_vs_i.ensembl-matched.degenes.fdr.05

rsem-prepare-reference-global-asm:

	cd assembly/global_merged; cat transcripts.fa.clean.nr | python $(protocol)/prepare-transcripts.py transcripts.fa.clean.nr.rsem knownIsoforms.txt
	cd assembly/global_merged; qsub -v "input=transcripts.fa.clean.nr.rsem,knownIsoforms=knownIsoforms.txt, \
		output=transcripts-rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-expression-global-asm:

	cd assembly/global_merged; qsub -v input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem",index="transcripts-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem",index="transcripts-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem",index="transcripts-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem",index="transcripts-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd assembly/global_merged; qsub -v input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem",index="transcripts-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem",index="transcripts-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="transcripts-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="transcripts-rsem" $(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-global-asm:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results  \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results  \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd assembly/global_merged; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd assembly/global_merged; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-assembly:

	cd assembly/global_merged; \
		python $(protocol)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

run-blast-assembly-gallus:

	cd assembly/global_merged; \
		python ~/rnaseq-comp-protocol/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd assembly/global_merged; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" ~/rnaseq-comp-protocol/blast.sh

run-blast-assembly-human:

	mkdir assembly/global_merged/Human_blast

	cd assembly/global_merged; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" ~/rnaseq-comp-protocol/blast.sh

	cd assembly/global_merged; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" ~/rnaseq-comp-protocol/blast.sh

run-goseq-assembly-gallus:

	Rcript $(protocol)/goseq_assembly_gallus.R

#############################
##### Cufflinks de novo #####
#############################

run-tophat-pe:

	cd tophat; qsub -v left=../reads/line6u.pe.1,right=../reads/line6u.pe.2,outdir=line6u_pe,index=gal4selected \
		$(protocol)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line6i.pe.1,right=../reads/line6i.pe.2,outdir=line6i_pe,index=gal4selected \
		$(protocol)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line7u.pe.1,right=../reads/line7u.pe.2,outdir=line7u_pe,index=gal4selected \
		$(protocol)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line7i.pe.1,right=../reads/line7i.pe.2,outdir=line7i_pe,index=gal4selected \
		$(protocol)/tophat_pe_job.sh

run-tophat-se:

	cd tophat; qsub -v input=../reads/line6u.se.fq,outdir=line6u_se,index=gal4selected \
		$(protocol)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line6i.se.fq,outdir=line6i_se,index=gal4selected \
		$(protocol)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line7u.se.fq,outdir=line7u_se,index=gal4selected \
		$(protocol)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line7i.se.fq,outdir=line7i_se,index=gal4selected \
		$(protocol)/tophat_se_job.sh

run-cufflinks:

	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		$(protocol)/cufflinks_job.sh; echo $$d; done

run-cuffmerge:

	cd tophat; cuffmerge -o merged_cuff_denovo -s gal4selected.fa -p 4 merge_list.txt

run-rsem-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; cat merged.gtf | python $(protocol)col/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_denovo; ~/rsem-1.2.7/rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-denovo
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-line7-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_denovo; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd tophat/merged_cuff_denovo; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
		python $(protocol)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 merged-denovo.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-cufflinks-denovo-gallus:

	cd tophat/merged_cuff_denovo; \
		python $(protocol)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
		python $(protocol)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd tophat/merged_cuff_denovo; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol)/blast.sh
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)/blast.sh

run-blast-cufflinks-denovo-human:

	mkdir tophat/merged_cuff_denovo/Human_blast
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol)/blast.sh
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)/blast.sh

####################################
###### Cufflinks + Ensembl #########
####################################

run-cuffmerge-ref:

	cd tophat; cuffmerge -o merged_cuffref --ref-gtf $(protocol)/Gallus_gallus.Galgal4.73.removed.gtf -s ../galGal4-removed.fa -p 4 merge_list.txt

run-rsem-cufflinks-ref:

	cd tophat/merged_cuff_ref; cat merged.gtf | python $(protocol)/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_ref; rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-ref
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \

		qsub -v index="merged-ref",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		$(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-line7-cufflinks-ref:

	cd tophat/merged_cuff_ref; \
	rsem-generate-data-matrix line7u-single-rsem-full.genes.results \
		line7u-paired-rsem-full.genes.results line7i-single-rsem-full.genes.results \
		line7i-paired-rsem-full.genes.results > line7u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_ref; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd tophat/merged_cuff_ref; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-cuffref:

	cd tophat/merged_cuff_ref; \
		python ~/rnaseq-comp-protocol/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 merged-ref.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

run-blast-cuffref-gallus:

	cd tophat/merged_cuff_ref; \
		python ~/rnaseq-comp-protocol/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd tophat/merged_cuff_ref; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" ~/rnaseq-comp-protocol/blast.sh

run-blast-cuffref-human:

	mkdir tophat/merged_cuff_ref/Human_blast
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" ~/rnaseq-comp-protocol/blast.sh

annotate-cuffref:

	cd tophat/merged_cuff_ref; \
		python $(protocol)/gene-rep-cufflinks.py merged-ref.transcripts.fa > merged-ref.genes.fa
	cd tophat/merged_cuff_ref; \
		qsub -v db="Gallus_prot",input="merged-ref.genes.fa",program="blastx",output="merged-ref.genes.xml" \
		$(protocol)/blast.sh

rsem-prepare-reference-cuffref-ensembl-matched:

	# cd tophat/merged_cuff_ref; \
	# 	python $(protocol)/get_top_hits.py merged-ref.genes.xml > merged-ref.tophits.txt
	# cd tophat/merged_cuff_ref; \
	# 	cat merged.rsem.gtf | python $(protocol)/cufflinks-to-known-isoforms.py > knownIsoforms.txt

	cd tophat/merged_cuff_ref; \
		python $(protocol)/get_best_ensembl_hits_cufflinks.py merged-ref.tophits.txt merged-ref.genes.fa \
		knownIsoforms.txt > merged.ensembl-matched.fa

	cd tophat/merged_cuff_ref; qsub -v \
	"input=merged.ensembl-matched.fa,knownIsoforms=knownIsoforms.txt,output=merged.ensembl-matched-rsem" \
	$(protocol)/rsem_prepare_reference.sh

run-rsem-cufflinks-ref-ensembl-matched:

	cd tophat/merged_cuff_ref; \
		qsub -v index="merged.ensembl-matched-rsem",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-ensembl-matched-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged.ensembl-matched-rsem",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-ensembl-matched-rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \
		qsub -v index="merged.ensembl-matched-rsem",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged.ensembl-matched-rsem",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-ensembl-matched-rsem" $(protocol)/rsem_calculate_expr_paired.sh

run-ebseq-cufflinks-ref-ensembl-matched:

	cd tophat/merged_cuff_ref; \
	rsem-generate-data-matrix line7u-single-ensembl-matched-rsem.genes.results \
		line7u-paired-ensembl-matched-rsem.genes.results \
		line7i-single-ensembl-matched-rsem.genes.results \
		line7i-paired-ensembl-matched-rsem.genes.results > line7u_vs_i.gene-ensembl-matched.counts.matrix
	cd tophat/merged_cuff_ref; rsem-run-ebseq line7u_vs_i.gene-ensembl-matched.counts.matrix 2,2 \
		line7u_vs_i.ensembl-matched.degenes
	cd tophat/merged_cuff_ref; rsem-control-fdr line7u_vs_i.ensembl-matched.degenes 0.05 \
		line7u_vs_i.ensembl-matched.degenes.fdr.05

run-goseq-cuffref-gallus:

	Rscript $(protocol)/goseq_cufflinks_gallus.R

########################################
######     Combined models    ##########
########################################

local-assembly: run-tophat-pe run-tophat-se extract-reads
run-tophat-pe:

	cd tophat; qsub -v outdir="line6u_pe",index="gal4selected",left="../reads/line6u.1_trim1.fastq",right="../reads/line6u.1_trim2.fastq",unpaired="../reads/line6u.1_trim_unpaired.fastq" $(protocol)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line6i_pe",index="gal4selected",left="../reads/line6i.1_trim1.fastq",right="../reads/line6i.1_trim2.fastq",unpaired="../reads/line6i.1_trim_unpaired.fastq" $(protocol)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7u_pe",index="gal4selected",left="../reads/line7u.1_trim1.fastq",right="../reads/line7u.1_trim2.fastq",unpaired="../reads/line7u.1_trim_unpaired.fastq" $(protocol)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7i_pe",index="gal4selected",left="../reads/line7i.1_trim1.fastq",right="../reads/line7i.1_trim2.fastq",unpaired="../reads/line7i.1_trim_unpaired.fastq" $(protocol)/tophat_pe_job.sh 

run-tophat-se:

	cd tophat; qsub -v outdir="line6u_se",index="gal4selected",input="../reads/line6u.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line6i_se",index="gal4selected",input="../reads/line6i.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7u_se",index="gal4selected",input="../reads/line7u.fq_trim.fastq" $(protocol)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7i_se",index="gal4selected",input="../reads/line7i.fq_trim.fastq" $(protocol)/tophat_se_job.sh

extract-reads:

	cd tophat; for dir in line??_?e; \
		do $(protocol)/extract_reads.sh $$dir/accepted_hits.bam $(protocol)/chromosomes.txt; \
	done

merge-bams:

	cd tophat; \
	for chr in $(cat $(protocol)/chromosomes.txt); do printf "merging %s..\n" "$chr";  \
		samtools merge -n merged/"$chr".bam \
		line6u_pe/"$chr".bam line6u_se/"$chr".bam \
		line6i_pe/"$chr".bam line6i_se/"$chr".bam \
		line7u_pe/"$chr".bam line7u_se/"$chr".bam \
		line7i_pe/"$chr".bam line7i_se/"$chr".bam; \
	done

run-velveth-local:

	cd tophat/merged; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" $(protocol)/velveth_local_job.sh; \
	done

run-velvetg-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol)/velvetg_local_job.sh; \
	done

run-oases-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol)/oases_local_job.sh; \
	done

combine-transcripts:

	cd tophat/merged; \
		for d in chr*_[0-9][0-9]; \
			do python ~/gimme/src/utils/rename_fasta.py $$d/transcripts.fa local_$$d >> local_merged.fa; \
	done

	cd assembly; \
		for d in global_[0-9][0-9]; \
			do python ~/gimme/src/utils/rename_fasta.py $$d/transcripts.fa global_$$d >> global_merged.fa; \
	done

clean-transcripts:

	cd tophat/merged; ~/seqclean-x86_64/seqclean local_merged.fa
	cd assembly; ~/seqclean-x86_64/seqclean global_merged.fa

remove-redundant-seq:

	cat tophat/merged/local_merged.fa.clean assembly/global_merged.fa.clean >> all.fa.clean
	qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh

	cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh

align-transcripts:

	python $(protocol)/split-fa.py all.fa.clean.nr
	for f in subsets*.fa; do \
		qsub -v input="$$f" $(protocol)/blat_job.sh; \
	done
	cat subsets*.fa.psl > all.fa.clean.nr.psl
	sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
	rm subsets*.fa.psl
	rm subsets*.fa

# run-cufflinks:
# 
# 	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
# 		$(protocol)/cufflinks_job.sh; echo $$d; done

run-cuffmerge-ref:

	cd tophat; cuffmerge -g ../Gallus_UCSC_ensembl_73.gtf.removed -o merged_cuff_ref -s gal4selected.fa -p 4 $(protocol)/merge_list.txt

# build-gene-models:
# 
# 	qsub -v input="all.fa.clean.nr.psl.best",ref="tophat/gal4selected.fa" $(protocol)/run_gimme.sh

# build-gene-models-with-cufflinks:
# 
# 	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_denovo/transcripts.gtf > tophat/merged_cuff_denovo/transcripts.bed
# 	qsub -v input1="all.fa.clean.nr.psl.best",input2="tophat/merged_cuff_denovo/transcripts.bed",ref="tophat/gal4selected.fa" $(protocol)/run_gimme2.sh

build-gene-models-with-cufflinks-ref:

	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_ref/merged.gtf > tophat/merged_cuff_ref/merged.bed
	cd combined; \
		qsub -v output="asm_cuff_ref_models.bed",input1="../all.fa.clean.nr.psl.best",input2="../tophat/merged_cuff_ref/merged.bed",ref="../tophat/gal4selected.fa" $(protocol)/run_gimme2.sh

models-to-transcripts:

	cd combined; python ~/gimme/src/utils/get_transcript_seq.py asm_cuff_ref_models.bed ../tophat/gal4selected.fa | sed 1d > asm_cuff_ref_models.bed.fa

rsem-combined-models-ref:

	cd combined; \
		cat asm_cuff_ref_models.bed.fa | python $(protocol)/fasta-to-gene-list.py > asm_cuff_ref_models.txt
	cd combined; \
		qsub -v list="asm_cuff_ref_models.txt",input="asm_cuff_ref_models.bed.fa",sample="asm_cuff_ref_models_rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-gimme-models-ref:

	cd combined; qsub -v input_read="../reads/line6u.se.fq",sample_name="line6u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line6i.se.fq",sample_name="line6i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol)/rsem_calculate_expr_single.sh

	cd combined; qsub -v input_read1="../reads/line6u.pe.1",input_read2="../reads/line6u.pe.2",sample_name="line6u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line6i.pe.1",input_read2="../reads/line6i.pe.2",sample_name="line6i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line7u.pe.1",input_read2="../reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line7i.pe.1",input_read2="../reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol)/rsem_calculate_expr_paired.sh

ebseq-line7-models-ref:

	cd combined; rsem-generate-data-matrix line7u-single-rsem-cuffref.genes.results \
		line7u-paired-rsem-cuffref.genes.results line7i-single-rsem-cuffref.genes.results \
		line7i-paired-rsem-cuffref.genes.results > line7u_vs_i.gene.cuffref.counts.matrix
	cd combined; rsem-run-ebseq line7u_vs_i.gene.cuffref.counts.matrix 2,2 line7u_vs_i.cuffref.degenes
	cd combined; rsem-control-fdr line7u_vs_i.cuffref.degenes 0.05 line7u_vs_i.cuffref.degenes.fdr.05

get-DE-sequences-combined:

	cd combined; \
		python $(protocol)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-combined:

	cd combined; \
		estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-combined-gallus:

	cd combined; python $(protocol)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
	cd combined; python $(protocol)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd combined; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol)/blast.sh
	cd combined; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)/blast.sh

run-blast-combined-human:

	mkdir assembly/global_merged/Human_blast
	cd combined; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol)/blast.sh
	cd combined; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol)/blast.sh


annotate-gimme:

	# python $(protocol)/get_top_hits.py asm_cuff_ref_models.longest.xml > asm_cuff_ref_models.longest.tophits.txt
	python $(protocol)/get_best_ensembl_hits_combined.py asm_cuff_ref_models.longest.tophits.txt \
		asm_cuff_ref_models.bed.fa > asm_cuff_ref_models.ensembl-matched.fa

	cat asm_cuff_ref_models.ensembl-matched.fa | python $(protocol)/fasta-to-gene-list.py \
		> asm_cuff_ref_models.ensembl-matched.txt
	qsub -v "knownIsoforms=asm_cuff_ref_models.ensembl-matched.txt,input=asm_cuff_ref_models.ensembl-matched.fa,\
		output=asm_cuff_ref_models_ensembl_matched_rsem" $(protocol)/rsem_prepare_reference.sh

rsem-calc-gimme-ensembl-matched:

	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref-ensembl-matched",index="asm_cuff_ref_models_ensembl_matched_rsem" $(protocol)/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref-ensembl-matched",index="asm_cuff_ref_models_ensembl_matched_rsem" $(protocol)/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref-ensembl-matched",index="asm_cuff_ref_models_ensembl_matched_rsem" $(protocol)/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref-ensembl-matched",index="asm_cuff_ref_models_ensembl_matched_rsem" $(protocol)/rsem_calculate_expr_paired.sh

ebseq-gimme-ensembl-matched:

	rsem-generate-data-matrix line7u-single-rsem-cuffref-ensembl-matched.genes.results \
		line7u-paired-rsem-cuffref-ensembl-matched.genes.results \
		line7i-single-rsem-cuffref-ensembl-matched.genes.results \
		line7i-paired-rsem-cuffref-ensembl-matched.genes.results > \
		line7u_vs_i.gene.cuffref-ensembl-matched.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.cuffref-ensembl-matched.counts.matrix 2,2 \
		line7u_vs_i.cuffref-ensembl-matched.degenes
	rsem-control-fdr line7u_vs_i.cuffref-ensembl-matched.degenes 0.05 \
		line7u_vs_i.cuffref-ensembl-matched.degenes.fdr.05

run-goseq-ensembl-human:

	cd human/ensembl; \
		Rscript $(protocol)/goseq_ensembl_human.R

run-goseq-ensembl-gallus:

	cd gallus/ensembl; \
	Rscript $(protocol)/goseq_ensembl_gallus.R

run-goseq-cufflinks-gallus:

	cd gallus/cufflinks; \
	Rscript $(protocol)/goseq_cufflinks_gallus.R

run-goseq-assembly-gallus:

	cd gallus/assembly; \
	Rscript $(protocol)/goseq_assembly_gallus.R

run-goseq-combined-gallus:

	cd gallus/combined; \
	Rscript $(protocol)/goseq_combined_gallus.R

create-combined-annotation:

	cd human_gallus; \
	Rscript $(protocol)/combine_kegg_annots.R; \
	python $(protocol)/create_kegg_annots.py line7u_vs_i.degenes.cuffref.tophits.nucl.gg.annots.txt \
		line7u_vs_i.degenes.cuffref.tophits.nucl.gg.pathways.txt line7u_vs_i.degenes.cuffref.tophits.nucl.hs.annots.txt \
		line7u_vs_i.degenes.cuffref.tophits.nucl.hs.pathways.txt ensembl_gallus_kegg_pathways.txt > \
		line7u_vs_i.degenes.cuffref.tophits.nucl.hs.gg.pathways.txt

run-goseq-combined-custom-gallus:

	Rscript $(protocol)/goseq_combined_custom_pathways.R
