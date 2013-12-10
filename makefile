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
