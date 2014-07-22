library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.ensembl.degenes.fdr.05',
                          stringsAsFactors=F, sep="\t", header=T)

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%rownames(degenes.table))
names(gene.vector)<-allgenes

categories <- read.table('custom_annotations.txt',
                         colClasses=c("character"),
                         sep='\t', header=FALSE)
colnames(categories) = c("geneID", "PATH")

pwf=nullp(gene.vector, 'galGal4', 'ensGene')

kegg = goseq(pwf, "galGal4", "ensGene",
             test.cats="KEGG", gene2cat=categories)

# Adjust P-value using BH method
kegg$padjust = p.adjust(kegg$over_represented_pvalue, method="BH")

kegg.sig = kegg[kegg$padjust<0.05,]

write.table(kegg.sig, 'line7u_vs_i.ensembl.degenes.custom.kegg.txt',
            sep='\t', row.names=F, quote=F)
