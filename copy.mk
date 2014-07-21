copypath=preeyano@hpc.msu.edu:/mnt/ls12/preeyanon/rnaseqcomp

ensembl:

	scp $(copypath)/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.ensembl.degenes.fdr.05

assembly:

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.assembly.degenes.fdr.05.tophits.gga

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.assembly.degenes.fdr.05.tophits.hsa

	scp $(copypath)/assembly/global_merged/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.assembly.degenes.fdr.05

cufflinks:

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.cufflinks.degenes.fdr.05.tophits.gga

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.cufflinks.degenes.fdr.05.tophits.hsa

	scp $(copypath)/tophat/merged_cuff_ref/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.cufflinks.degenes.fdr.05

gimme:

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.gga.tophits \
		results/line7u_vs_i.gimme.degenes.fdr.05.tophits.gga

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05.hsa.tophits \
		results/line7u_vs_i.gimme.degenes.fdr.05.tophits.hsa

	scp $(copypath)/gimme/line7u_vs_i.degenes.fdr.05 \
		results/line7u_vs_i.gimme.degenes.fdr.05
