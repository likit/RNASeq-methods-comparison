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

run-blast-human-full:

	# for a custom pathway annotation
	python $(protocol)/gene-rep-ensbl.py Gallus_gallus.Galgal4.73.cdna.all.fa \
		> ensembl.fa.longest

	qsub -v "db=Human_prot,input=ensembl.fa.longest,program=blastx,\
		output=ensembl.fa.longest.hsa.xml" $(protocol)/blast.sh

run-blast-human:

	# create a directory for human annotation
	qsub -v "db=Human_prot,input=line7u_vs_i.degenes.fdr.05.fa.longest,\
		program=blastx,output=human/line7u_vs_i.degenes.fdr.05.fa.longest.hsa.xml" \
		$(protocol)/blast.sh

get-top-hits-human:

	python $(protocol)/get_top_hits.py ensembl.fa.longest.hsa.xml \
	> all-ensembl-hsa-tophits.txt

get-tophits-degenes:

	python $(protocol)/get_top_hits.py line7u_vs_i.degenes.fdr.05.fa.longest.hsa.xml \
	> ensembl-hsa-tophits.txt

	python $(protocol)/tophits-to-degenes-ensembl.py \
		line7u_vs_i.degenes.fdr.05 ensembl-hsa-tophits.txt > \
		line7u_vs_i.ensembl.degenes.fdr.05.tophits.hsa
