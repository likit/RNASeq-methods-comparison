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

All software for Linux 64-bit machine can be downloaded at
http://athyra.ged.msu.edu/~preeyano/software/.

#Data processing

Set up protocol and gimme path:

    export PROTOCOL=<path to protocol root directory>
    export GIMMEDIR=<path to Gimme root directory>

###Cufflinks

Run RSEM prepare reference:

    make -f $PROTOCOL/cufflinks.mk prepare-reference-cufflinks

Run RSEM:

    make -f $PROTOCOL/cufflinks.mk run-rsem-calc-expression-cufflinks

Run EBseq:

    make -f $PROTOCOL/cufflinks.mk run-ebseq-cufflinks

###Global assembly


Run RSEM prepare reference:

    make -f $PROTOCOL/global_assembly.mk rsem-prepare-reference-assembly

Run RSEM:

    make -f $PROTOCOL/global_assembly.mk rsem-calc-expression-assembly

Run EBseq:

    make -f $PROTOCOL/global_assembly.mk run-ebseq-assembly

###Merged models

Build merged models:

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL gimmedir=$GIMMEDIR build-merged-gene-models

Note, gimme.bed may contain warning messages from pygr at the beginning of the file.
The messages have to be removed before running the next step.

Build RSEM reference for Ensembl-matched merged models:

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL rsem-prepare-reference-merged-models

Run RSEM on Ensembl-matched merged models:

    make -f $PROTOCOL/gimme.mk protocol=$PROTOCOL rsem-calc-expression-merged-models

Run EBseq:

    make -f $PROTOCOL/gimme.mk ebseq-gimme
