get-all-pathways:

	cd results; Rscript $(protocol)/get_all_ensembl_pathways.R

run-goseq-ensembl-gallus:

	cd results; Rscript $(protocol)/goseq_ensembl_gallus.R

run-goseq-gimme-gallus:

	cd results; Rscript $(protocol)/goseq_combined_gallus.R

run-goseq-assembly-human:

	# cd human/assembly; \
	# 	python $(protocol)/get_top_hits.py \
	# 	line7u_vs_i.degenes.fdr.05.fa.longest.xml > \
	# 	line7u_vs_i.degenes.tophits.nucl.txt

	cd human/assembly; \
		Rscript $(protocol)/goseq_assembly_human.R

run-goseq-cufflinks-human:

	cd human/cufflinks; \
		python $(protocol)/get_top_hits.py \
		line7u_vs_i.cuffref.degenes.fdr.05.fa.longest.xml > \
		line7u_vs_i.degenes.tophits.nucl.txt

	cd human/cufflinks; \
		Rscript $(protocol)/goseq_assembly_human.R

run-goseq-combined-human:

	cd human/combined; \
		python $(protocol)/get_top_hits.py \
		line7u_vs_i.cuffref.degenes.fdr.05.fa.longest.xml > \
		line7u_vs_i.degenes.tophits.nucl.txt

	cd human/combined; \
		Rscript $(protocol)/goseq_assembly_human.R
