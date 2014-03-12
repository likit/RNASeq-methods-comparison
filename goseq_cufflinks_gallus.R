library(goseq)
library(org.Gg.eg.db)
library(KEGG.db)
library(ggplot2)
library(biomaRt)

degenes.table<-read.table('line7u_vs_i.cuffref.tophits.nucl.txt', stringsAsFactors=F, sep="\t", header=T)
annots<-select(org.Gg.eg.db, keys=degenes.table$geneID, columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")

annotated.degenes<-merge(degenes.table, annots, by.x="geneID", by.y="ENSEMBL")

# remove duplicated Entrez ID
uniq.annotated.degenes<-annotated.degenes[!duplicated(annotated.degenes$geneID),]

# remove gene with no Entrez ID
uniq.annotated.degenes<-uniq.annotated.degenes[!is.na(uniq.annotated.degenes$ENTREZID),]

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
allgenes<-getBM(attributes='ensembl_gene_id', mart=mart)
allgenes<-allgenes$ensembl_gene_id


gene.vector<-as.integer(allgenes%in%degenes.table$geneID)
names(gene.vector)<-allgenes

pwf=nullp(gene.vector, 'galGal4', 'ensGene')

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
get_genes_kegg = function(cat, data, prefix)
{
    m = match(xx[[cat]], data$ENTREZID)
    mm = m[!is.na(m)]
    d = data.frame(cat, data[mm,]$geneID, data[mm,]$ENTREZID)

    filename = paste(prefix, cat, sep="_")
    write.table(d, filename, sep="\t", row.names=F, col.names=F, quote=F)
    return(d)
}
df = lapply(KEGG_SIG$category, get_genes_kegg, uniq.annotated.degenes, "line7_goseq_KEGG_genes")

# Get number of genes in each pathway
#KEGG_SIG$ngenes = sapply(df, nrow) # get a number of genes from each category

# Writing pathway information to a file
cat("Writing pathways' info to a file...\n")
write.table(KEGG_SIG, 'line7u_vs_i.cuffref.degenes.KEGG.txt', sep='\t', row.names=F, quote=F)

write.table(uniq.annotated.degenes, 'uniq-annotated-degenes.txt', sep='\t', row.names=F, col.names=F, quote=F)
# Cleveland plot
#cat("Plotting...\n")
#KEGG_SIG$log10padjust=(-1)*log10(KEGG_SIG$padjust)
#nameorder = KEGG_SIG$pathway[sort.list(KEGG_SIG$log10padjust, decreasing=F)]
#KEGG_SIG$pathway = factor(KEGG_SIG$pathway, levels=nameorder)
#ggplot(KEGG_SIG, aes(x=pathway, y=log10padjust, size=ngenes)) +
#    geom_point(colour="grey30") +
#    theme_bw() +
#    theme(panel.grid.major.y = element_blank(),
#        panel.grid.minor.y=element_blank(),
#        panel.grid.major.x=element_line(colour="grey60",
#        linetype="dashed"), axis.text.x=element_text(angle=90, hjust=1)) +
#    scale_size_area(max_size=14) +
#    labs(list(title="Line 7 Post Infection", y="log10(adjusted p-value)"))
#ggsave("line7_KEGG_cleveland_inverse.pdf")
#cat("Done!\n")
