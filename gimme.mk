# Requires Gimme
# Please set gimmetools to gimme/src/utils
#
combine-transcripts:

	cd local_assembly; \
		for d in chr*_[0-9][0-9]; \
			do python $(gimmetools)/rename_fasta.py $$d/transcripts.fa local_$$d >> local_merged.fa; \
	done

	cd assembly; \
		for d in global_[0-9][0-9]; \
			do python $(gimmetools)/rename_fasta.py $$d/transcripts.fa global_$$d >> global_merged.fa; \
	done

clean-transcripts:

	cd local_assembly; seqclean local_merged.fa -c 8
	cd assembly; seqclean global_merged.fa -c 8

remove-redundant-seq:

	if [ ! -d gimme ]; then mkdir gimme; fi
	cat local_assembly/local_merged.fa.clean assembly/global_merged.fa.clean >> gimme/all.fa.clean
	# cd gimme; \
	# 	qsub -v input="all.fa.clean",output="all.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh

	# cd assembly; qsub -v input="global_merged.fa.clean",output="global_merged.fa.clean.nr",c="1.0" $(protocol)/cdhit_job.sh

align-transcripts:

	cd gimme; \
		python $(protocol)/split-fa.py all.fa.clean.nr 10000 subsets
	cd gimme; \
		for f in subsets*.fa; do \
			qsub -v input="$$f" $(protocol)/blat_job.sh; \
		done
		# wait till all the jobs are done before proceeding

sort-alignments:

	cd gimme; \
		cat subsets*.fa.psl > all.fa.clean.nr.psl; \
		sort -k 10 all.fa.clean.nr.psl > all.fa.clean.nr.psl.sorted; \
		pslReps -nohead -singleHit all.fa.clean.nr.psl.sorted all.fa.clean.nr.psl.best info; \
		rm subsets*.fa.psl; \
		rm subsets*.fa

build-merged-gene-models:

	cd gimme; \
	python ~/gimme/src/utils/gff2bed.py ../tophat/merged_cuff_ref/merged.gtf > \
		../tophat/merged_cuff_ref/merged.bed; \
	qsub -v "output=gimme-models.bed,input1=all.fa.clean.nr.psl.best,\
		input2=../tophat/merged_cuff_ref/merged.bed,\
		ref=../gal4selected.fa" $(protocol)/run_gimme.sh
