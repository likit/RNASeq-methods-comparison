get-all-pathways:

	cd results; Rscript $(protocol)/get_all_ensembl_pathways.R

create-custom-annotations:

	# create chicken-human annotation
	cd results; \
	python $(protocol)/ensembl-human-tophits-to-annots.py \
		all-ensembl-hsa-tophits.txt > all-ensembl-hsa-annotations.txt

	cd results; \
	python $(protocol)/create_kegg_annots_ensembl.py \
		all-ensembl-hsa-annotations.txt ensembl_gallus_kegg_pathways.txt \
		ensembl_human_kegg_pathways.txt > custom_annotations.txt

run-goseq-ensembl-gallus:

	cd results; Rscript $(protocol)/goseq_ensembl_gallus.R

run-goseq-ensembl-human:

	cd results; Rscript $(protocol)/goseq_ensembl_human.R

run-goseq-ensembl-custom:

	cd results; Rscript $(protocol)/goseq_ensembl_custom.R

run-goseq-gimme-gallus:

	cd results; Rscript $(protocol)/goseq_combined_gallus.R

run-goseq-gimme-human:

	cd results; Rscript $(protocol)/goseq_combined_human.R

run-goseq-gimme-custom:

	cd results; Rscript $(protocol)/goseq_combined_custom.R

# run-goseq-cufflinks-gallus:
# 
# 	cd results; Rscript $(protocol)/goseq_cufflinks_gallus.R
# 
# run-goseq-cufflinks-human:
# 
# 	cd results; Rscript $(protocol)/goseq_cufflinks_human.R

run-goseq-cufflinks-custom:

	cd results; Rscript $(protocol)/goseq_cufflinks_custom.R

# run-goseq-assembly-gallus:
# 
# 	cd results; Rscript $(protocol)/goseq_assembly_gallus.R
# 
# run-goseq-assembly-human:
# 
# 	cd results; Rscript $(protocol)/goseq_assembly_human.R

run-goseq-assembly-custom:

	cd results; Rscript $(protocol)/goseq_assembly_custom.R
