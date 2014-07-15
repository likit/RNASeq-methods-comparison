##############################################
############# Ensembl annotations #############
###############################################

prepare-ensembl-reference:

	python $(protocol)/gtf_to_known_isoforms.py \
		Gallus_gallus.Galgal4.73.removed.gtf > Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt

	rsem-prepare-reference --gtf Gallus_gallus.Galgal4.73.removed.gtf \
		--transcript-to-gene-map Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt \
		gal4selected.fa ensembl_genes

calc-expression:

	qsub -v "input_read=reads/line7u.se.fq,sample_name=line7u-single-rsem,\
		index=ensembl_genes" $(protocol)/rsem_calculate_expr_single.sh
	qsub -v "input_read=reads/line7i.se.fq,sample_name=line7i-single-rsem,\
		index=ensembl_genes" $(protocol)/rsem_calculate_expr_single.sh

	qsub -v "input_read1=reads/line7u.pe.1,input_read2=reads/line7u.pe.2,\
		sample_name=line7u-paired-rsem,index=ensembl_genes" \
		$(protocol)/rsem_calculate_expr_paired.sh
	qsub -v "input_read1=reads/line7i.pe.1,input_read2=reads/line7i.pe.2,\
		sample_name=line7i-paired-rsem,index=ensembl_genes" \
		$(protocol)/rsem_calculate_expr_paired.sh
