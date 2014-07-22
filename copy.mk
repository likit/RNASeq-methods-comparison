copypath=preeyano@hpc.msu.edu:/mnt/ls12/preeyanon/rnaseqcomp

misc:

	# scp $(copypath)/ensembl-hsa-tophits.txt\
	# 	results/ensembl-hsa-tophits.txt

	# scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.gga.tophits \
	# 	results/line7u_vs_i.gimme.degenes.fdr.05.gga.tophits

	# scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.hsa.tophits \
	# 	results/line7u_vs_i.gimme.degenes.fdr.05.hsa.tophits

	scp $(copypath)/all-ensembl-hsa-tophits.txt \
		results/all-ensembl-hsa-tophits.txt

ensembl:

	scp $(copypath)/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.ensembl.degenes.fdr.05
	scp $(copypath)/line7u_vs_i.ensembl.degenes.fdr.05.tophits.hsa \
		results/line7u_vs_i.ensembl.degenes.fdr.05.tophits.hsa

	scp $(copypath)/line7u-single-rsem.genes.results \
		results/line7u-single-rsem.genes.results.ensembl
	scp $(copypath)/line7i-single-rsem.genes.results \
		results/line7i-single-rsem.genes.results.ensembl

	scp $(copypath)/line7u-paired-rsem.genes.results \
		results/line7u-paired-rsem.genes.results.ensembl
	scp $(copypath)/line7i-paired-rsem.genes.results \
		results/line7i-paired-rsem.genes.results.ensembl

assembly:

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.assembly.degenes.fdr.05.tophits.gga

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.assembly.degenes.fdr.05.tophits.hsa

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.assembly.degenes.fdr.05

	scp $(copypath)/assembly/global_merged/line7u-single-rsem.genes.results \
		results/line7u-single-rsem.genes.results.assembly
	scp $(copypath)/assembly/global_merged/line7i-single-rsem.genes.results \
		results/line7i-single-rsem.genes.results.assembly

	scp $(copypath)/assembly/global_merged/line7u-paired-rsem.genes.results \
		results/line7u-paired-rsem.genes.results.assembly
	scp $(copypath)/assembly/global_merged/line7i-paired-rsem.genes.results \
		results/line7i-paired-rsem.genes.results.assembly

cufflinks:

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.cufflinks.degenes.fdr.05.tophits.gga

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.cufflinks.degenes.fdr.05.tophits.hsa

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.cufflinks.degenes.fdr.05

	scp $(copypath)/tophat/merged_cuff_ref/line7u-single-rsem.genes.results \
		results/line7u-single-rsem.genes.results.cufflinks
	scp $(copypath)/tophat/merged_cuff_ref/line7i-single-rsem.genes.results \
		results/line7i-single-rsem.genes.results.cufflinks

	scp $(copypath)/tophat/merged_cuff_ref/line7u-paired-rsem.genes.results \
		results/line7u-paired-rsem.genes.results.cufflinks
	scp $(copypath)/tophat/merged_cuff_ref/line7i-paired-rsem.genes.results \
		results/line7i-paired-rsem.genes.results.cufflinks

gimme:

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.gimme.degenes.fdr.05.tophits.gga

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.gimme.degenes.fdr.05.tophits.hsa

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.gimme.degenes.fdr.05

	scp $(copypath)/gimme/line7u-single-rsem.genes.results \
		results/line7u-single-rsem.genes.results.gimme
	scp $(copypath)/gimme/line7i-single-rsem.genes.results \
		results/line7i-single-rsem.genes.results.gimme

	scp $(copypath)/gimme/line7u-paired-rsem.genes.results \
		results/line7u-paired-rsem.genes.results.gimme
	scp $(copypath)/gimme/line7i-paired-rsem.genes.results \
		results/line7i-paired-rsem.genes.results.gimme
