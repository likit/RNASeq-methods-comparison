library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table(
                  'gallus/combined/line7u_vs_i.cuffref.tophits.nucl.txt',
                  stringsAsFactors=F, sep="\t", header=T)

annots<-select(org.Gg.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots,
                         by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[
                          !duplicated(annotated.degenes$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[
                           !is.na(uniq.annotated.degenes$ENTREZID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")

allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%degenes.table$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')
categories = read.table(
    'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.hs.gg.pathways.txt')

colnames(categories) = c("geneID", "PATH")

KEGG = goseq(pwf, "galGal4", "ensGene",
             test.cats="KEGG", gene2cat=categories)

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.05,]
KEGG_SIG$category=sprintf("%05s", KEGG_SIG$category)

write.table(KEGG_SIG,
  'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.hs.gg.KEGG.txt',
  quote=F, row.names=F, sep='\t')
