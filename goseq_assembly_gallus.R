library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.assembly.degenes.fdr.05.tophits.gga',
                          stringsAsFactors=F, sep="\t", header=T)
colnames(degenes.table)<-c("seqID", "geneID")
annots<-select(org.Gg.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID","PATH"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots,
                         by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[
                          !duplicated(annotated.degenes$geneID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%degenes.table$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')

kegg = goseq(pwf, "galGal4", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
kegg$padjust = p.adjust(kegg$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
kegg.sig = kegg[kegg$padjust<0.05,]
pathway = stack(mget(kegg.sig$category, KEGGPATHID2NAME))
kegg.sig$pathway = pathway$values

write.table(kegg.sig, 'line7u_vs_i.assembly.degenes.gga.kegg.txt', sep='\t',
            row.names=F, quote=F)
write.table(uniq.annotated.degenes[!is.na(uniq.annotated.degenes$PATH),],
            'line7u_vs_i.assembly.degenes.gga.kegg.id.txt', sep='\t',
            row.names=F, col.names=F, quote=F)
