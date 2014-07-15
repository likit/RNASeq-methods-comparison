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

run-ebseq-line7:

	# requires RSEM 1.2.7
	rsem-generate-data-matrix line7u-single-rsem.genes.results \
		line7u-paired-rsem.genes.results line7i-single-rsem.genes.results \
		line7i-paired-rsem.genes.results > line7u_vs_i.gene.counts.matrix
	rsem-run-ebseq line7u_vs_i.gene.counts.matrix 2,2 line7u_vs_i.degenes
	rsem-control-fdr line7u_vs_i.degenes 0.05 line7u_vs_i.degenes.fdr.05

get-de-seqs:

	python $(protocol)/rsem-output-ensbl-to-fasta.py line7u_vs_i.degenes.fdr.05 \
		ensembl_genes.transcripts.fa Gallus_gallus.Galgal4.73.removed.knownIsoforms.txt \
		> line7u_vs_i.degenes.fdr.05.fa

run-blast-gallus:

	python $(protocol)/gene-rep-ensbl.py line7u_vs_i.degenes.fdr.05.fa \
		> line7u_vs_i.degenes.fdr.05.fa.longest

	qsub -v "db=Gallus_prot,input=line7u_vs_i.degenes.fdr.05.fa.longest, \
		program=blastx,output=line7u_vs_i.degenes.fdr.05.fa.longest.xml" \
		$(protocol)/blast.sh

run-blast-human:

	# create a directory for human annotation
	if [ ! -d human ]; then mkdir human;  fi
	qsub -v "db=Human_prot,input=line7u_vs_i.degenes.fdr.05.fa.longest, \
		program=blastx,output=human/line7u_vs_i.degenes.fdr.05.fa.longest.xml" \
		$(protocol)/blast.sh
