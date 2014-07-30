RNASeq-methods-comparison
=========================

This project aims to compare the results of biological network analysis
from RNASeq data using gene models from following methods:
* Tophat/Cufflinks (with and without integration of Ensembl gene models)
* _De novo_ assembly
* _De novo_ assembly + local assembly (reference guided)
* Ensembl predicted gene models

Required software
-----------------

Please make sure all required software is available on your machine.
* Bowtie 1.0 for RSEM
* Bowtie 2.1.0 for Tophat
* Tophat 2.0.9
* Cufflinks 2.1.1
* Velvet 1.2.03
* Oases 0.2.06
* BLAT
* RSEM 1.2.7
* Condetri 2.1 for quality trimming
* Seqclean
* Biopython

Required R packages
-------------------

All packages are available at bioconductor.org.

* org.Gg.eg.db
* org.Hs.eg.db
* KEGG.db
* GO.db
* biomaRt
* goseq

All software for Linux 64-bit machine can be downloaded at
http://athyra.ged.msu.edu/~preeyano/software/.

Set up protocol and gimme path

    export PROTOCOL=<path to protocol root directory>
    export GIMMEDIR=<path to Gimme root directory>

###Prepare Reads

    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-pe
    make -f $PROTOCOL/misc.mk protocol=$PROTOCOL run-quality-trim-se

###Cufflinks

Map reads to the genome using Tophat

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-tophat-pe
    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-tophat-se

Run RSEM prepare reference

    make -f $PROTOCOL/cufflinks.mk prepare-reference-cufflinks

Then run Cufflinks

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-cufflinks

Merge gene models from all samples

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL run-cuffmerge-ref

Annotate genes with chicken and human Ensembl models

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL annotate-cufflinks

Prepare RSEM reference

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL prepare-reference-cufflinks

Run RSEM

    make -f $PROTOCOL/cufflinks.mk run-rsem-calc-expression-cufflinks

Run EBseq with FDR 0.05

    make -f $PROTOCOL/cufflinks.mk run-ebseq-cufflinks

Get top hits for DE genes from chicken and human:

    make -f $PROTOCOL/cufflinks.mk protocol=$PROTOCOL get-tophits-degenes

###Global assembly

Run RSEM prepare reference

    make -f $PROTOCOL/global_assembly.mk rsem-prepare-reference-assembly

Run RSEM

    make -f $PROTOCOL/global_assembly.mk rsem-calc-expression-assembly

Run EBseq with FDR 0.05

    make -f $PROTOCOL/global_assembly.mk run-ebseq-assembly

Get top hits for DE genes from chicken and human

    make -f $PROTOCOL/global_assembly.mk protocol=$PROTOCOL get-tophits-degenes

###Local assembly

Extract reads from each chromosome

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL extract-reads

Merge reads from each chromosome to local assembly directory

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL merge-bams

Run Velveth, Velvetg and Oases

    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-velveth-local
    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-velvetg-local
    make -f $PROTOCOL/local_assembly.mk protocol=$PROTOCOL run-oases-local

Merge transcripts

    make -f $PROTOCOL/local_assembly.mk gimmedir=$GIMMEDIR combine-transcripts

Clean transcripts and remove redundant transcripts

    make -f $PROTOCOL/local_assembly.mk clean-transcripts
    make -f $PROTOCOL/local_assembly.mk remove-redundant-seq

Align all transcripts to the genome

    make -f $PROTOCOL/local_assembly.mk align-transcripts:

###Merged models

Build merged models

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL gimmedir=$GIMMEDIR build-merged-gene-models

Note, gimme.bed may contain warning messages from pygr at the beginning of the file.
The messages have to be removed before running the next step.

Build RSEM reference for Ensembl-matched merged models

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL rsem-prepare-reference-merged-models

Run RSEM on Ensembl-matched merged models

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL rsem-calc-expression-merged-models

Run EBseq with FDR 0.05

    make -f $PROTOCOL/gimme.mk ebseq-gimme

Get top hits for DE genes from chicken and human
