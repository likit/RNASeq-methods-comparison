extract-reads:

	cd local_assembly; for dir in ../tophat/line??_?e; \
		do $(protocol)/extract_reads.sh $$dir/accepted_hits.bam $(protocol)/chromosomes.txt; \
	done

merge-bams:

	if [ ! -d local_assembly ]; then mkdir local_assembly; fi
	cd tophat; \
	for chr in $$(cat $(protocol)/chromosomes.txt); do printf "merging %s..\n" "$$chr";  \
		samtools merge -n ../local_assembly/"$$chr".bam \
		line6u_pe/"$$chr".bam line6u_se/"$$chr".bam \
		line6i_pe/"$$chr".bam line6i_se/"$$chr".bam \
		line7u_pe/"$$chr".bam line7u_se/"$$chr".bam \
		line7i_pe/"$$chr".bam line7i_se/"$$chr".bam; \
	done
	cd local_assembly; rm chr32.bam # no reads mapped to chromosome 32

run-velveth-local:

	cd local_assembly; \
	for f in *.bam; \
		do qsub -v outdir=$$(basename "$$f" .bam),input="$$f" $(protocol)/velveth_local_job.sh; \
	done

run-velvetg-local:

	cd local_assembly; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol)/velvetg_local_job.sh; \
	done

run-oases-local:

	cd local_assembly; \
	for d in chr*_*; \
		do qsub -v indir="$$d" $(protocol)/oases_local_job.sh; \
	done

combine-transcripts:

	cd tophat/merged; \
	for d in chr*_[0-9][0-9]; \
		do python $(gimmedir)/src/utils/rename_fasta.py $$d/transcripts.fa local_$$d >> local_merged.fa; \
	done

	cd assembly; \
	for d in global_[0-9][0-9]; \
		do python $(gimmedir)/src/utils/rename_fasta.py $$d/transcripts.fa global_$$d >> global_merged.fa; \
	done

clean-transcripts:

	cd local_assembly; seqclean local_merged.fa
	cd assembly; seqclean global_merged.fa

remove-redundant-seq:

	cat tophat/merged/local_merged.fa.clean assembly/global_merged.fa.clean >> all.fa.clean
	qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh

	cd assembly; qsub -v "input=global_merged.fa.clean,\
		output=global_merged.fa.clean.nr,c=1.0" $(protocol)/cdhit_job.sh

align-transcripts:

	python $(protocol)/split-fa.py all.fa.clean.nr

	qsub -v "input=all.fa.clean.nr" $(protocol)/blat_job.sh; done

	sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted
	pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info
