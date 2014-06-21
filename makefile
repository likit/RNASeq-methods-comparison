protocol_path = ~/projects/rnaseq-compare/protocol

###############################################
############# Ensembl annotations #############
###############################################

preprocess-ucsc-gene-models:

	cat Gallus_UCSC_ensembl_73.gtf | cut -f 1 | sort | uniq -c | \
		grep -v -e random -e chrUn > Gallus_UCSC_ensembl_73.gtf.removed
	cat Gallus_UCSC_ensembl_73.txt | cut -f 1 | sort | uniq -c | \
		grep -v -e random -e chrUn > Gallus_UCSC_ensembl_73.txt.removed
	cat Gallus_UCSC_ensembl_73.txt.removed | cut -f 2,13 | \
		grep -v name | awk -v OFS="\t" '{print $2,$1}' > Gallus_UCSC_ensembl_73.knownIsoforms.txt

run-rsem-prepare-reference:

	rsem-prepare-reference --gtf Gallus_UCSC_ensembl_73.gtf.removed \
		--transcript-to-gene-map Gallus_UCSC_ensembl_73.knownIsoforms.txt \
		galGal4-removed.fa galGal4-removed

run-rsem-calc-expression:

	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh

run-ebseq-line7:

	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-ensembl:

	python $(protocol_path)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 merged-ref.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-ensembl:

	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol_path)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-ensembl-gallus:

	python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
	python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
	qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-blast-ensembl-human:

	mkdir Human_blast
	qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
	qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-goseq-ensembl-gallus:

	Rscript $(protocol_path)/goseq_ensembl_gallus.R

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

	cd assembly; qsub -v pe_input="paired.fastq",se_input="single.fastq" $(protocol_path)/velveth_job.sh

run-velvetg:

	cd assembly; qsub $(protocol_path)/velvetg_job.sh

run-oases:

	cd assembly; qsub $(protocol_path)/oases_job.sh

run-oasesM:

	cd assembly; qsub $(protocol_path)/velvethM_job.sh
	cd assembly; qsub $(protocol_path)/velvetgM_job.sh
	cd assembly; qsub $(protocol_path)/oasesM_job.sh

clean-transcripts:

	# -A needed to keep poly-A tail
	cd assembly/global_merged; ~/seqclean-x86_64/seqclean transcripts.fa -c 8 -A -o transcripts.fa.clean
	qsub -v input="assembly/global_merged/transcripts.fa.clean",output="assembly/global_merged/transcripts.fa.clean.nr",c="1.0" protocol/cdhit_job.sh

run-rsem-prepare-reference-global-asm:

	cd assembly/global_merged; cat transcripts.fa.clean.nr | python $(protocol_path)/prepare-transcripts.py transcripts.fa.clean.nr.rsem knownIsoforms.txt
	cd assembly/global_merged; qsub $(protocol_path)/rsem_prepare_reference.sh

run-rsem-calc-expression-global-asm:

	cd assembly/global_merged; qsub -v input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem",index="transcripts-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem",index="transcripts-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem",index="transcripts-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem",index="transcripts-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh

	cd assembly/global_merged; qsub -v input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem",index="transcripts-rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem",index="transcripts-rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="transcripts-rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="transcripts-rsem" $(protocol_path)/rsem_calculate_expr_paired.sh

run-ebseq-line7-global-asm:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results  \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results  \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd assembly/global_merged; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd assembly/global_merged; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-assembly:

	cd assembly/global_merged; \
		python $(protocol_path)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-assembly:

	cd assembly/global_merged; \
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol_path)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-assembly-gallus:

	cd assembly/global_merged; \
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd assembly/global_merged; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-blast-assembly-human:

	mkdir assembly/global_merged/Human_blast
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-goseq-assembly-gallus:

	Rcript $(protocol_path)/goseq_assembly_gallus.R

#############################
##### Cufflinks de novo #####
#############################

run-tophat-pe:

	cd tophat; qsub -v left=../reads/line6u.pe.1,right=../reads/line6u.pe.2,outdir=line6u_pe,index=gal4selected \
		$(protocol_path)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line6i.pe.1,right=../reads/line6i.pe.2,outdir=line6i_pe,index=gal4selected \
		$(protocol_path)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line7u.pe.1,right=../reads/line7u.pe.2,outdir=line7u_pe,index=gal4selected \
		$(protocol_path)/tophat_pe_job.sh
	cd tophat; qsub -v left=../reads/line7i.pe.1,right=../reads/line7i.pe.2,outdir=line7i_pe,index=gal4selected \
		$(protocol_path)/tophat_pe_job.sh

run-tophat-se:

	cd tophat; qsub -v input=../reads/line6u.se.fq,outdir=line6u_se,index=gal4selected \
		$(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line6i.se.fq,outdir=line6i_se,index=gal4selected \
		$(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line7u.se.fq,outdir=line7u_se,index=gal4selected \
		$(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v input=../reads/line7i.se.fq,outdir=line7i_se,index=gal4selected \
		$(protocol_path)/tophat_se_job.sh

run-cufflinks:

	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		$(protocol_path)/cufflinks_job.sh; echo $$d; done

run-cuffmerge:

	cd tophat; cuffmerge -o merged_cuff_denovo -s gal4selected.fa -p 4 merge_list.txt

run-rsem-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; cat merged.gtf | python $(protocol_path)col/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_denovo; ~/rsem-1.2.7/rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-denovo
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh

run-ebseq-line7-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_denovo; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd tophat/merged_cuff_denovo; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
		python $(protocol_path)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 merged-denovo.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol_path)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-cufflinks-denovo-gallus:

	cd tophat/merged_cuff_denovo; \
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd tophat/merged_cuff_denovo; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-blast-cufflinks-denovo-human:

	mkdir tophat/merged_cuff_denovo/Human_blast
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

####################################
###### Cufflinks + Ensembl #########
####################################

run-cuffmerge-ref:

	cd tophat; cuffmerge -o merged_cuff_ref --ref-gtf ../Gallus_UCSC_ensembl_73.gtf -s gal4selected.fa -p 4 merge_list.txt

run-rsem-cufflinks-ref:

	cd tophat/merged_cuff_ref; cat merged.gtf | python $(protocol_path)/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_ref; rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-ref
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \

		qsub -v index="merged-ref",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		$(protocol_path)/rsem_calculate_expr_paired.sh

run-ebseq-line7-cufflinks-ref:

	cd tophat/merged_cuff_ref; \
	rsem-generate-data-matrix line7u-single-rsem-full.genes.results \
		line7u-paired-rsem-full.genes.results line7i-single-rsem-full.genes.results \
		line7i-paired-rsem-full.genes.results > line7u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_ref; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd tophat/merged_cuff_ref; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-DE-sequences-cuffref:

	cd tophat/merged_cuff_ref; \
		python $(protocol_path)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 merged-ref.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-cuffref:

	cd tophat/merged_cuff_ref; \
	estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol_path)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-cuffref-gallus:

	cd tophat/merged_cuff_ref; \
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
		python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd tophat/merged_cuff_ref; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-blast-cuffref-human:

	mkdir tophat/merged_cuff_ref/Human_blast
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-goseq-cuffref-gallus:

	Rscript $(protocol_path)/goseq_cufflinks_gallus.R

########################################
######     Combined models    ##########
########################################

local-assembly: run-tophat-pe run-tophat-se extract-reads
run-tophat-pe:

	cd tophat; qsub -v outdir="line6u_pe",index="gal4selected",left="../reads/line6u.1_trim1.fastq",right="../reads/line6u.1_trim2.fastq",unpaired="../reads/line6u.1_trim_unpaired.fastq" $(protocol_path)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line6i_pe",index="gal4selected",left="../reads/line6i.1_trim1.fastq",right="../reads/line6i.1_trim2.fastq",unpaired="../reads/line6i.1_trim_unpaired.fastq" $(protocol_path)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7u_pe",index="gal4selected",left="../reads/line7u.1_trim1.fastq",right="../reads/line7u.1_trim2.fastq",unpaired="../reads/line7u.1_trim_unpaired.fastq" $(protocol_path)/tophat_pe_job.sh 
	cd tophat; qsub -v outdir="line7i_pe",index="gal4selected",left="../reads/line7i.1_trim1.fastq",right="../reads/line7i.1_trim2.fastq",unpaired="../reads/line7i.1_trim_unpaired.fastq" $(protocol_path)/tophat_pe_job.sh 

run-tophat-se:

	cd tophat; qsub -v outdir="line6u_se",index="gal4selected",input="../reads/line6u.fq_trim.fastq" $(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line6i_se",index="gal4selected",input="../reads/line6i.fq_trim.fastq" $(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7u_se",index="gal4selected",input="../reads/line7u.fq_trim.fastq" $(protocol_path)/tophat_se_job.sh
	cd tophat; qsub -v outdir="line7i_se",index="gal4selected",input="../reads/line7i.fq_trim.fastq" $(protocol_path)/tophat_se_job.sh

extract-reads:

	cd tophat; for dir in line??_?e; \
		do $(protocol_path)/extract_reads.sh $$dir/accepted_hits.bam $(protocol_path)/chromosomes.txt; \
	done

merge-bams:

	cd tophat; \
	for chr in $(cat $(protocol_path)/chromosomes.txt); do printf "merging %s..\n" "$chr";  \
		samtools merge -n merged/"$chr".bam \
		line6u_pe/"$chr".bam line6u_se/"$chr".bam \
		line6i_pe/"$chr".bam line6i_se/"$chr".bam \
		line7u_pe/"$chr".bam line7u_se/"$chr".bam \
		line7i_pe/"$chr".bam line7i_se/"$chr".bam; \
	done

run-velveth-local:

	cd tophat/merged; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" $(protocol_path)/velveth_local_job.sh; \
	done

run-velvetg-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol_path)/velvetg_local_job.sh; \
	done

run-oases-local:

	cd tophat/merged; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol_path)/oases_local_job.sh; \
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
	qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" $(protocol_path)/cdhit_job.sh

	cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" $(protocol_path)/cdhit_job.sh

align-transcripts:

	python $(protocol_path)/split-fa.py all.fa.clean.nr
	for f in subsets*.fa; do \
		qsub -v input="$$f" $(protocol_path)/blat_job.sh; \
	done
	cat subsets*.fa.psl > all.fa.clean.nr.psl
	sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
	rm subsets*.fa.psl
	rm subsets*.fa

# run-cufflinks:
# 
# 	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
# 		$(protocol_path)/cufflinks_job.sh; echo $$d; done

run-cuffmerge-ref:

	cd tophat; cuffmerge -g ../Gallus_UCSC_ensembl_73.gtf.removed -o merged_cuff_ref -s gal4selected.fa -p 4 $(protocol_path)/merge_list.txt

# build-gene-models:
# 
# 	qsub -v input="all.fa.clean.nr.psl.best",ref="tophat/gal4selected.fa" $(protocol_path)/run_gimme.sh

# build-gene-models-with-cufflinks:
# 
# 	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_denovo/transcripts.gtf > tophat/merged_cuff_denovo/transcripts.bed
# 	qsub -v input1="all.fa.clean.nr.psl.best",input2="tophat/merged_cuff_denovo/transcripts.bed",ref="tophat/gal4selected.fa" $(protocol_path)/run_gimme2.sh

build-gene-models-with-cufflinks-ref:

	python ~/gimme/src/utils/gff2bed.py tophat/merged_cuff_ref/merged.gtf > tophat/merged_cuff_ref/merged.bed
	cd combined; \
		qsub -v output="asm_cuff_ref_models.bed",input1="../all.fa.clean.nr.psl.best",input2="../tophat/merged_cuff_ref/merged.bed",ref="../tophat/gal4selected.fa" $(protocol_path)/run_gimme2.sh

models-to-transcripts:

	cd combined; python ~/gimme/src/utils/get_transcript_seq.py asm_cuff_ref_models.bed ../tophat/gal4selected.fa | sed 1d > asm_cuff_ref_models.bed.fa

rsem-combined-models-ref:

	cd combined; \
		cat asm_cuff_ref_models.bed.fa | python $(protocol_path)/fasta-to-gene-list.py > asm_cuff_ref_models.txt
	cd combined; \
		qsub -v list="asm_cuff_ref_models.txt",input="asm_cuff_ref_models.bed.fa",sample="asm_cuff_ref_models_rsem" $(protocol_path)/rsem_prepare_reference.sh

rsem-calc-gimme-models-ref:

	cd combined; qsub -v input_read="../reads/line6u.se.fq",sample_name="line6u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line6i.se.fq",sample_name="line6i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line7u.se.fq",sample_name="line7u-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh
	cd combined; qsub -v input_read="../reads/line7i.se.fq",sample_name="line7i-single-rsem-cuffref",index="asm_cuff_ref_models_rsem" \
		$(protocol_path)/rsem_calculate_expr_single.sh

	cd combined; qsub -v input_read1="../reads/line6u.pe.1",input_read2="../reads/line6u.pe.2",sample_name="line6u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line6i.pe.1",input_read2="../reads/line6i.pe.2",sample_name="line6i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line7u.pe.1",input_read2="../reads/line7u.pe.2",sample_name="line7u-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol_path)/rsem_calculate_expr_paired.sh
	cd combined; qsub -v input_read1="../reads/line7i.pe.1",input_read2="../reads/line7i.pe.2",sample_name="line7i-paired-rsem-cuffref",index="asm_cuff_ref_models_rsem" $(protocol_path)/rsem_calculate_expr_paired.sh

ebseq-line7-models-ref:

	cd combined; rsem-generate-data-matrix line7u-single-rsem-cuffref.genes.results \
		line7u-paired-rsem-cuffref.genes.results line7i-single-rsem-cuffref.genes.results \
		line7i-paired-rsem-cuffref.genes.results > line7u_vs_i.gene.cuffref.counts.matrix
	cd combined; rsem-run-ebseq line7u_vs_i.gene.cuffref.counts.matrix 2,2 line7u_vs_i.cuffref.degenes
	cd combined; rsem-control-fdr line7u_vs_i.cuffref.degenes 0.05 line7u_vs_i.cuffref.degenes.fdr.05

get-DE-sequences-combined:

	cd combined; \
		python $(protocol_path)/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa

translate-DE-sequences-combined:

	cd combined; \
		estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M $(protocol_path)/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

run-blast-combined-gallus:

	cd combined; python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest
	cd combined; python $(protocol_path)/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.longest

	cd combined; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
	cd combined; \
		qsub -v db="Gallus_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-blast-combined-human:

	mkdir assembly/global_merged/Human_blast
	cd combined; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml" $(protocol_path)/blast.sh
	cd combined; \
		qsub -v db="Human_prot",input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="Human_blast/line7u_vs_i.degenes.fdr.05.fa.longest.xml" $(protocol_path)/blast.sh

run-goseq-ensembl-gallus:

	cd gallus/ensembl; \
	Rscript $(protocol_path)/goseq_ensembly_gallus.R

run-goseq-cufflinks-gallus:

	cd gallus/cufflinks; \
	Rscript $(protocol_path)/goseq_cufflinks_gallus.R

run-goseq-assembly-gallus:

	cd gallus/assembly; \
	Rscript $(protocol_path)/goseq_assembly_gallus.R

run-goseq-combined-gallus:

	cd gallus/combined; \
	Rscript $(protocol_path)/goseq_combined_gallus.R

create-combined-annotation:

	cd human_gallus; \
	Rscript $(protocol_path)/combine_kegg_annots.R; \
	python $(protocol_path)/create_kegg_annots.py line7u_vs_i.degenes.cuffref.tophits.nucl.gg.annots.txt \
		line7u_vs_i.degenes.cuffref.tophits.nucl.gg.pathways.txt line7u_vs_i.degenes.cuffref.tophits.nucl.hs.annots.txt \
		line7u_vs_i.degenes.cuffref.tophits.nucl.hs.pathways.txt ensembl_gallus_kegg_pathways.txt > \
		line7u_vs_i.degenes.cuffref.tophits.nucl.hs.gg.pathways.txt

run-goseq-combined-custom-gallus:

	Rscript $(protocol_path)/goseq_combined_custom_pathways.R
