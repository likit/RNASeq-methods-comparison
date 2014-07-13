library("org.Gg.eg.db")
library("org.Hs.eg.db")

degenes.gg.table<-read.table('gallus/gimme/line7u_vs_i.cuffref.tophits.nucl.txt',
                          stringsAsFactors=F, sep="\t", header=T)
ggannots<-select(org.Gg.eg.db,
               keys=degenes.table$geneID,
               columns=c("SYMBOL","ENTREZID", "PATH"),
               keytype="ENSEMBL")
degenes.hs.table<-read.table('human/gimme/line7u_vs_i.cuffref.tophits.nucl.txt',
                          stringsAsFactors=F, sep="\t", header=T)
hsannots<-select(org.Hs.eg.db,
                 keys=degenes.table$geneID,
                 columns=c("SYMBOL","ENTREZID", "PATH"),
                 keytype="ENSEMBL")
write.table(ggannots,
            'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.gg.annots.txt',
            sep='\t', quote=F, row.names=F, col.names=F)
write.table(hsannots,
            'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.hs.annots.txt',
            sep='\t', quote=F, row.names=F, col.names=F)
write.table(data.frame(rownames(degenes.gg.table), degenes.gg.table$geneID),
            'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.gg.txt',
            sep='\t', quote=F, row.names=F, col.names=F)
write.table(data.frame(rownames(degenes.hs.table), degenes.hs.table$geneID),
            'human_gallus/line7u_vs_i.degenes.cuffref.tophits.nucl.hs.txt',
            sep='\t', quote=F, row.names=F, col.names=F)