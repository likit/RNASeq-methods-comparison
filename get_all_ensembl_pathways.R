library("biomaRt")
library("org.Gg.eg.db")
library("org.Hs.eg.db")
library("KEGG.db")

# Get all chicken KEGG pathways
mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
all.gg.genes<-getBM(attributes='ensembl_gene_id', mart=mart)
all.gg.genes<-all.gg.genes$ensembl_gene_id

gg.annots<-select(org.Gg.eg.db, keys=all.gg.genes,
               columns=c("SYMBOL","ENTREZID", "PATH"), keytype="ENSEMBL")
write.table(gg.annots, 'ensembl_gallus_kegg_pathways.txt',
            sep='\t', col.names=F, row.names=F, quote=F)

gg.pathway = stack(mget(gg.annots[!is.na(gg.annots$PATH),]$PATH,
                        KEGGPATHID2NAME))
write.table(gg.pathway, 'all_gallus_pathway_names.txt',
            sep='\t', quote=F, row.names=F, col.names=F)

# Get all human KEGG pathways
mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
all.hs.genes<-getBM(attributes='ensembl_gene_id', mart=mart)
all.hs.genes<-all.hs.genes$ensembl_gene_id

hs.annots<-select(org.Hs.eg.db, keys=all.hs.genes,
                  columns=c("SYMBOL","ENTREZID", "PATH"),
                  keytype="ENSEMBL")
write.table(hs.annots, 'ensembl_human_kegg_pathways.txt',
            sep='\t', col.names=F, row.names=F, quote=F)

hs.pathway = stack(mget(hs.annots[!is.na(hs.annots$PATH),]$PATH,
                        KEGGPATHID2NAME))
write.table(hs.pathway, 'all_human_pathway_names.txt',
            sep='\t', quote=F, row.names=F, col.names=F)
