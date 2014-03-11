##### RSEM from Ensembl annotations #####
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

	qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="galGal4-removed" \
		protocol/rsem_calculate_expr_paired.sh

ebseq-line6:

	rsem-generate-data-matrix line6u-single-rsem.genes.results \
		line6u-paired-rsem.genes.results line6i-single-rsem.genes.results \
		line6i-paired-rsem.genes.results > line6u_vs_i.gene.counts.matrix
	rsem-run-ebseq line6u_vs_i.gene.counts.matrix 2,2 line6u_vs_i.degenes
	rsem-control-fdr line6u_vs_i.degenes 0.05 line6u_vs_i.degenes.fdr.05

ebseq-line7:

	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

prepare-for-blast2go:

	#python protocol/rsem-output-ensbl-to-fasta.py line6u_vs_i.degenes.fdr.05 Gallus_gallus.Galgal4.73.pep.all.fa > line6u_vs_i.degenes.fdr.05.fa 
	#python protocol/rsem-output-ensbl-to-fasta.py line7u_vs_i.degenes.fdr.05 Gallus_gallus.Galgal4.73.pep.all.fa > line7u_vs_i.degenes.fdr.05.fa 
	#python protocol/gene-rep-ensbl.py line6u_vs_i.degenes.fdr.05.fa > line6u_vs_i.degenes.fdr.05.fa.prot.longest
	#python protocol/gene-rep-ensbl.py line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.prot.longest
	#python protocol/split-fa.py line6u_vs_i.degenes.fdr.05.fa.prot.longest 100 line6u_vs_i.degenes.fdr.05.fa.prot.longest
	#python protocol/split-fa.py line7u_vs_i.degenes.fdr.05.fa.prot.longest 100 line7u_vs_i.degenes.fdr.05.fa.prot.longest
	
	mkdir Human_blast
	#for f in *prot.longest*.fa; do \
	#	qsub -v input="$$f",program="blastp" ~/rnaseq-comp-protocol/blast.sh; \
	#done
	qsub -v input="line6u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp" ~/rnaseq-comp-protocol/blast.sh; \
	qsub -v input="line7u_vs_i.degenes.fdr.05.fa.prot.longest",program="blastp" ~/rnaseq-comp-protocol/blast.sh; \

run-blast2go-human:

	qsub -v input="Human_blast/line7u_vs_i.degenes.fdr.05.fa.prot.longest.xml",outdir="Human_blast/line7u_blast2go_outdir" ~/rnaseq-comp-protocol/b2g_job.sh
	qsub -v input="Human_blast/line6u_vs_i.degenes.fdr.05.fa.prot.longest.xml",outdir="Human_blast/line6u_blast2go_outdir" ~/rnaseq-comp-protocol/b2g_job.sh

run-interpro-ensembl:

	#mkdir interpro
	for f in *prot*longest*fa; do \
		qsub -v input="$$f" ~/rnaseq-comp-protocol/iprscan.sh; \
	done

run-blast2go:

	qsub -v input="line7u_vs_i.degenes.fdr.05.fa.prot.longest.fa.xml",outdir="line7u_blast2go_outdir" protocol/b2g_job.sh
	qsub -v input="line6u_vs_i.degenes.fdr.05.fa.prot.longest.fa.xml",outdir="line6u_blast2go_outdir" protocol/b2g_job.sh

##### De novo assembly (global + local) #####

run-quality-trim-pe:

	perl ~/condetri_v2.1.pl -fastq1=reads/line7u.pe.1 -fastq2=reads/line7u.pe.2 -cutfirst 10 -sc=33
	qsub -v left=reads/line7i.pe.1,right=reads/line7i.pe.2 protocol/quality_trim_pe_job.sh

run-quality-trim-se:

	for r in reads/*.se.fq; do qsub -v input="$$r" protocol/quality_trim_se_job.sh; done

interleave-reads:

	cd assembly; ~/velvet_1.2.03/shuffleSequences_fastq.pl pe.1.fastq pe.2.fastq paired.fastq

run-velveth:

	cd assembly; qsub -v pe_input="paired.fastq",se_input="single.fastq" ../protocol/velveth_job.sh

run-velvetg:

	cd assembly; qsub ../protocol/velvetg_job.sh

run-oases:

	cd assembly; qsub ../protocol/oases_job.sh

run-oasesM:

	cd assembly; qsub ../protocol/velvethM_job.sh
	cd assembly; qsub ../protocol/velvetgM_job.sh
	cd assembly; qsub ../protocol/oasesM_job.sh

#De novo assembly only
######################

clean-transcripts:
	# -A needed to keep poly-A tail
	cd assembly/global_merged; ~/seqclean-x86_64/seqclean transcripts.fa -c 8 -A -o transcripts.fa.clean
	qsub -v input="assembly/global_merged/transcripts.fa.clean",output="assembly/global_merged/transcripts.fa.clean.nr",c="1.0" protocol/cdhit_job.sh

run-rsem-prepare-reference-global-asm:

	cd assembly/global_merged; cat transcripts.fa.clean.nr | python ../../protocol/prepare-transcripts.py transcripts.fa.clean.nr.rsem knownIsoforms.txt
	cd assembly/global_merged; qsub ../../protocol/rsem_prepare_reference.sh

run-rsem-calc-expression-global-asm:

	cd assembly/global_merged; qsub -v input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem",index="transcripts-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem",index="transcripts-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem",index="transcripts-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd assembly/global_merged; qsub -v input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem",index="transcripts-rsem" \
		../../protocol/rsem_calculate_expr_single.sh

	cd assembly/global_merged; qsub -v input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem",index="transcripts-rsem" ../../protocol/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem",index="transcripts-rsem" ../../protocol/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem",index="transcripts-rsem" ../../protocol/rsem_calculate_expr_paired.sh
	cd assembly/global_merged; qsub -v input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem",index="transcripts-rsem" ../../protocol/rsem_calculate_expr_paired.sh

ebseq-line6-global-asm:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line6u-single-rsem.genes.results \
		line6u-paired-rsem.genes.results line6i-single-rsem.genes.results \
		line6i-paired-rsem.genes.results > line6u_vs_i.gene.counts.matrix
	cd assembly/global_merged; rsem-run-ebseq line6u_vs_i.gene.counts.matrix 2,2 line6u_vs_i.degenes
	cd assembly/global_merged; rsem-control-fdr line6u_vs_i.degenes 0.05 line6u_vs_i.degenes.fdr.05

ebseq-line7-global-asm:

	cd assembly/global_merged; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results  \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results  \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd assembly/global_merged; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd assembly/global_merged; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

prepare-for-blast2go-global-asm:

	#cd assembly/global_merged; \
	#estscan -t line6u_vs_i.degenes.fdr.05.fa.prot -M ../../protocol/gallus.hm line6u_vs_i.degenes.fdr.05.fa > line6u_vs_i.degenes.fdr.05.fa.nucl

	#cd assembly/global_merged; \
	#estscan -t line7u_vs_i.degenes.fdr.05.fa.prot -M ../../protocol/gallus.hm line7u_vs_i.degenes.fdr.05.fa > line7u_vs_i.degenes.fdr.05.fa.nucl

	#cd assembly/global_merged; \
	#	python ../../protocol/rsem-output-to-fasta.py line7u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line7u_vs_i.degenes.fdr.05.fa
	#cd assembly/global_merged; \
	#	python ../../protocol/rsem-output-to-fasta.py line6u_vs_i.degenes.fdr.05 transcripts-rsem.transcripts.fa > line6u_vs_i.degenes.fdr.05.fa
	#cd assembly/global_merged; \
	#	python ../../protocol/gene-rep-velvet.py line6u_vs_i.degenes.fdr.05.fa.prot > line6u_vs_i.degenes.fdr.05.fa.prot.longest
	#cd assembly/global_merged; \
	#	python ../../protocol/gene-rep-velvet.py line7u_vs_i.degenes.fdr.05.fa.prot > line7u_vs_i.degenes.fdr.05.fa.prot.longest

	#cd assembly/global_merged; \
	#	python ../../protocol/split-fa.py line7u_vs_i.degenes.fdr.05.fa.prot.longest 100 line7u_vs_i.degenes.fdr.05.fa.prot.longest
	#cd assembly/global_merged; \
	#	python ../../protocol/split-fa.py line6u_vs_i.degenes.fdr.05.fa.prot.longest 100 line6u_vs_i.degenes.fdr.05.fa.prot.longest

	#mkdir assembly/global_merged/Human_blast
	#cd assembly/global_merged; \
	#for f in *prot.longest*.fa; do \
	#	qsub -v input="$$f",program="blastp" ../../protocol/blast.sh; \
	#done

	cd assembly/global_merged; qsub -v input="line6u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="human_blast" ~/rnaseq-comp-protocol/blast.sh
	cd assembly/global_merged; qsub -v input="line7u_vs_i.degenes.fdr.05.fa.longest",program="blastx",output="human_blast" ~/rnaseq-comp-protocol/blast.sh

run-interpro-global-asm:

	cd assembly/global_merged/interpro; \
		for f in ../*prot*longest*fa; do \
			qsub -v input="$$f" ../../../protocol/iprscan.sh; \
		done

##### Tuxedo suit (Tophat + Cufflinks) #####

run-tophat-pe:
	cd tophat; qsub -v left=../reads/line6u.pe.1,right=../reads/line6u.pe.2,outdir=line6u_pe,index=gal4selected \
		../protocol/tophat_pe_job.sh

run-tophat-se:

	cd tophat; qsub -v input=../reads/line7u.se.fq,outdir=line7u_se,index=gal4selected \
		../protocol/tophat_se_job.sh

run-cufflinks:

	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		../protocol/cufflinks_job.sh; echo $$d; done

run-cuffmerge:

	cd tophat; cuffmerge -o merged_cuff_denovo -s gal4selected.fa -p 4 merge_list.txt

run-rsem-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; cat merged.gtf | python ../../protocol/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_denovo; ~/rsem-1.2.7/rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-denovo
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; \
		qsub -v index="merged-denovo",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh

ebseq-line6-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	rsem-generate-data-matrix line6u-single-rsem.genes.results \
		line6u-paired-rsem.genes.results line6i-single-rsem.genes.results \
		line6i-paired-rsem.genes.results > line6u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_denovo; rsem-run-ebseq line6u_vs_i.gene.counts.matrix 2,2 line6u_vs_i.degenes
	cd tophat/merged_cuff_denovo; rsem-control-fdr line6u_vs_i.degenes 0.05 line6u_vs_i.degenes.fdr.05

ebseq-line7-cufflinks-denovo:

	cd tophat/merged_cuff_denovo; \
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	cd tophat/merged_cuff_denovo; rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	cd tophat/merged_cuff_denovo; rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

run-rsem-cufflinks-ref:

	cd tophat/merged_cuff_ref; cat merged.gtf | python ../../protocol/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_ref; ~/rsem-1.2.7/rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-ref
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		../../protocol/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_ref; \

		qsub -v index="merged-ref",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_ref; \
		qsub -v index="merged-ref",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		../../protocol/rsem_calculate_expr_paired.sh

