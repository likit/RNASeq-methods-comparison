library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.gimme.degenes.fdr.05.tophits.gga',
                          stringsAsFactors=F, sep="\t", header=T)

colnames(degenes.table)<-c("SeqId", "geneID")
annots<-select(org.Gg.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots, by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[
                          !duplicated(annotated.degenes$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[
                           !is.na(uniq.annotated.degenes$ENTREZID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")

allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%uniq.annotated.degenes$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')

KEGG = goseq(pwf, "galGal4", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.05,]
pathway = stack(mget(KEGG[KEGG$padjust<0.05,]$category, KEGGPATHID2NAME))
KEGG_SIG$pathway = pathway$values

write.table(KEGG_SIG, 'line7_vs_i.gimme.degenes.KEGG.txt',
            sep='\t', row.names=F, quote=F)

write.table(uniq.annotated.degenes, 'gimme.uniq-annotated-degenes.txt',
            sep='\t', row.names=F, col.names=F, quote=F)
