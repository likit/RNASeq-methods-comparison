library(goseq)
library(org.Hs.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.gimme.degenes.fdr.05.tophits.hsa',
                          stringsAsFactors=F, sep="\t", header=T)

colnames(degenes.table)<-c("seqId", "geneID")
annots<-select(org.Hs.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID", "PATH"),
               keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots,
                         by.x="geneID", by.y="ENSEMBL")

# remove duplicated genes
uniq.annotated.degenes<-annotated.degenes[
                          !duplicated(annotated.degenes$geneID),]

allgenes<-read.table('all-ensembl-hsa-annotations.txt',
                     sep='\t', header=F, stringsAsFactor=F)

allgenes<-sort(allgenes$V2)
allgenes<-allgenes[!duplicated(allgenes)]

gene.vector<-as.integer(allgenes%in%degenes.table$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'hg19', 'ensGene')

kegg = goseq(pwf, 'hg19', 'ensGene', test.cats="KEGG")

# Adjust P-value using BH method
kegg$padjust = p.adjust(kegg$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
kegg.sig = kegg[kegg$padjust<0.05,]
pathway = stack(mget(kegg.sig$category, KEGGPATHID2NAME))
kegg.sig$pathway = pathway$values

write.table(kegg.sig, 'line7u_vs_i.gimme.degenes.hsa.kegg.txt',
            sep='\t', row.names=F, quote=F)

write.table(uniq.annotated.degenes[!is.na(uniq.annotated.degenes$PATH),],
            'line7u_vs_i.gimme.degenes.hsa.kegg.id.txt',
            sep='\t', row.names=F, col.names=F, quote=F)
