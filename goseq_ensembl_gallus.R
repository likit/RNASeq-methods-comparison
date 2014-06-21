library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(ggplot2)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.degenes.fdr.05', stringsAsFactors=F, sep="\t", header=T)
annots<-select(org.Gg.eg.db, keys=rownames(degenes.table), columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots, by.x=0, by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[!duplicated(annotated.degenes$Row.names),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[!is.na(uniq.annotated.degenes$ENTREZID),]
mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%uniq.annotated.degenes$Row.names)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')
#GO.wall <- goseq(pwf, "hg19", "ensGene")
#enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]

cat("KEGG pathway analysis..\n")
# KEGG Pathway analysis
KEGG = goseq(pwf, "galGal4", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.05,]
pathway = stack(mget(KEGG[KEGG$padjust<0.05,]$category, KEGGPATHID2NAME))
KEGG_SIG$pathway = pathway$values
xx = as.list(org.Gg.egPATH2EG)
xx = xx[!is.na(xx)] # remove KEGG IDs that do not match any gene

cat("Writing genes to files..\n")
# Write genes in each pathway to separate files
get_genes_kegg = function(c, data, prefix)
{
    m = match(xx[[c]], data$ENTREZID)
    mm = m[!is.na(m)]
    d = data.frame(c, data[mm,]$Row.names, data[mm,]$ENTREZID)

    filename = paste(prefix, c, sep="_")
    write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
    return(d)
}
df = lapply(KEGG_SIG$category, get_genes_kegg, uniq.annotated.degenes, "line7_goseq_KEGG_genes")

# Get number of genes in each pathwayKEGG_SIG$ngenes = sapply(df, nrow) # get a number of genes from each category

# Writing pathway information to a file
cat("Writing pathways' info to a file...\n")
write.table(KEGG_SIG, 'line7u_vs_i.degenes.KEGG.txt', sep='\t', row.names=F, quote=F)
write.table(uniq.annotated.degenes, 'uniq-annotated-degenes.txt', sep='\t', row.names=F, col.names=F, quote=F)

cat("Done!\n")
