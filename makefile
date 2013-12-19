preprocess-ucsc-gene-models:
	cat Gallus_UCSC_ensembl_73.gtf | cut -f 1 | sort | uniq -c | \
		grep -v -e random -e chrUn > Gallus_UCSC_ensembl_73.gtf.removed
	cat Gallus_UCSC_ensembl_73.txt | cut -f 1 | sort | uniq -c | \
		grep -v -e random -e chrUn > Gallus_UCSC_ensembl_73.txt.removed
	cat Gallus_UCSC_ensembl_73.txt.removed | cut -f 2,13 | \
		grep -v name | awk -v OFS="\t" '{print $2,$1}' > Gallus_UCSC_ensembl_73.knownIsoforms.txt

prepare-reference-rsem:
	rsem-prepare-reference --gtf Gallus_UCSC_ensembl_73.gtf.removed \
		--transcript-to-gene-map Gallus_UCSC_ensembl_73.knownIsoforms.txt \
		galGal4-removed.fa galGal4-removed

rsem-calc-expression:
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
