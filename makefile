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
	qsub -v input_read="reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		protocols/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		protocols/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		protocols/rsem_calculate_expr_single.sh
	qsub -v input_read="reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		protocols/rsem_calculate_expr_single.sh

	qsub -v input_read1="reads/line6u.pe.1",input_read2="reads/line6u.pe.2",sample_name="line6u-paired-rsem" \
		protocols/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line6i.pe.1",input_read2="reads/line6i.pe.2",sample_name="line6i-paired-rsem" \
		protocols/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7u.pe.1",input_read2="reads/line7u.pe.2",sample_name="line7u-paired-rsem" \
		protocols/rsem_calculate_expr_paired.sh
	qsub -v input_read1="reads/line7i.pe.1",input_read2="reads/line7i.pe.2",sample_name="line7i-paired-rsem" \
		protocols/rsem_calculate_expr_paired.sh

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

run-quality-trim-pe:
	# perl ~/condetri_v2.1.pl -fastq1=reads/line7u.pe.1 -fastq2=reads/line7u.pe.2 -cutfirst 10 -sc=33
	qsub -v left=reads/line7i.pe.1,right=reads/line7i.pe.2 protocols/quality_trim_pe_job.sh

run-quality-trim-se:
	for r in reads/*.se.fq; do qsub -v input="$$r" protocols/quality_trim_se_job.sh; done

interleave-reads:
	cd assembly; ~/velvet_1.2.03/shuffleSequences_fastq.pl pe.1.fastq pe.2.fastq paired.fastq

run-velveth:
	cd assembly; qsub -v pe_input="paired.fastq",se_input="single.fastq" ../protocols/velveth_job.sh

run-velvetg:
	cd assembly; qsub ../protocols/velvetg_job.sh

run-tophat-pe:
	cd tophat; qsub -v left=../reads/line6u.pe.1,right=../reads/line6u.pe.2,outdir=line6u_pe,index=gal4selected \
		../protocols/tophat_pe_job.sh

run-tophat-se:
	cd tophat; qsub -v input=../reads/line7u.se.fq,outdir=line7u_se,index=gal4selected \
		../protocols/tophat_se_job.sh

run-cufflinks:
	cd tophat; for d in line??_?e; do qsub -v outdir="$$d",input="$$d/accepted_hits.bam" \
		../protocols/cufflinks_job.sh; echo $$d; done

run-cuffmerge:
	cd tophat; cuffmerge -o merged_cuff_denovo -s gal4selected.fa -p 4 merge_list.txt

run-rsem-cufflinks-denovo:
	cd tophat/merged_cuff_denovo; cat merged.gtf | python ../../protocols/fix-gtf.py > merged.rsem.gtf 
	cd tophat/merged_cuff_denovo; ~/rsem-1.2.7/rsem-prepare-reference --gtf merged.rsem.gtf ../../galGal4-removed.fa merged-denovo
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read="../../reads/line6u.se.fq",sample_name="line6u-single-rsem" \
		../../protocols/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read="../../reads/line6i.se.fq",sample_name="line6i-single-rsem" \
		../../protocols/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read="../../reads/line7u.se.fq",sample_name="line7u-single-rsem" \
		../../protocols/rsem_calculate_expr_single.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read="../../reads/line7i.se.fq",sample_name="line7i-single-rsem" \
		../../protocols/rsem_calculate_expr_single.sh

	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read1="../../reads/line6u.pe.1",input_read2="../../reads/line6u.pe.2",sample_name="line6u-paired-rsem" ../../protocols/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read1="../../reads/line6i.pe.1",input_read2="../../reads/line6i.pe.2",sample_name="line6i-paired-rsem" ../../protocols/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read1="../../reads/line7u.pe.1",input_read2="../../reads/line7u.pe.2",sample_name="line7u-paired-rsem" ../../protocols/rsem_calculate_expr_paired.sh
	cd tophat/merged_cuff_denovo; qsub -v index="merged-denovo",input_read1="../../reads/line7i.pe.1",input_read2="../../reads/line7i.pe.2",sample_name="line7i-paired-rsem" ../../protocols/rsem_calculate_expr_paired.sh
