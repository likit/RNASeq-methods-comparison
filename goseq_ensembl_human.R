library(goseq)
library(org.Hs.eg.db)
library(KEGG.db)
library(ggplot2)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.degenes.tophits.prot.txt',
                          stringsAsFactors=F, sep="\t", header=T)
annots<-select(org.Hs.eg.db, keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots, by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[
                          !duplicated(annotated.degenes$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[
                       !is.na(uniq.annotated.degenes$ENTREZID),]

mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id

gene.vector<-as.integer(allgenes%in%uniq.annotated.degenes$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'hg19', 'ensGene')

# KEGG Pathway analysis
KEGG = goseq(pwf, "hg19", "ensGene", test.cats="KEGG")

# Adjust P-value using BH method
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue, method="BH")

# Get pathway names for significant patways
KEGG_SIG = KEGG[KEGG$padjust<0.05,]
pathway = stack(mget(KEGG[KEGG$padjust<0.05,]$category, KEGGPATHID2NAME))
KEGG_SIG$pathway = pathway$values
xx = as.list(org.Hs.egPATH2EG)
xx = xx[!is.na(xx)] # remove KEGG IDs that do not match any gene

cat("Writing genes to files..\n")
# Write genes in each pathway to separate files
get_genes_kegg = function(cat, data, prefix)
{
    m = match(xx[[cat]], data$ENTREZID)
    mm = m[!is.na(m)]
    d = data.frame(cat, data[mm,]$geneID, data[mm,]$ENTREZID)

    filename = paste(prefix, cat, sep="_")
    write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
    return(d)
}
df = lapply(KEGG_SIG$category, get_genes_kegg,
                uniq.annotated.degenes,
                "line7_goseq_KEGG_genes")

write.table(KEGG_SIG, 'line7u_vs_i.degenes.KEGG.txt', sep='\t',
            row.names=F, quote=F)
write.table(uniq.annotated.degenes, 'line7_vs_i.cuffref.degenes.with.kegg.id.txt',
            sep='\t', row.names=F, col.names=F, quote=F)
