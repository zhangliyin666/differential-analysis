install.packages('GOplot')
library(GOplot)
data(EC)
head(EC$david)
head(EC$genelist)


deg <- read.table(file = "sig_gene.txt", sep = "\t", header = T, stringsAsFactors = F)
which(is.na(deg),arr.ind = T)
deg[is.na(deg)] <- 0


deg <- deg[,c(8,3)]
colnames(deg) <- c("ID","logFC")

gene<-deg$ID

gene.df<-bitr(gene, fromType = "SYMBOL",
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Mm.eg.db)

ego_cc<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_bp<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf<-enrichGO(gene       = gene.df$SYMBOL,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)

ego_bp <- simplify(ego_bp,cutoff=0.7,by="p.adjust",select_fun=min)
ego_cc <- simplify(ego_cc,cutoff=0.7,by="p.adjust",select_fun=min)
ego_mf <- simplify(ego_mf,cutoff=0.7,by="p.adjust",select_fun=min)

ego1 <-data.frame(ego_bp)
ego2 <-data.frame(ego_cc)
ego3 <-data.frame(ego_mf)

ego1$Category <- "BP"
ego2$Category <- "CC"
ego3$Category <- "MF"

ego <- rbind(ego1,ego2,ego3)
ego <- ego[,c(10,1,2,8,6)]
colnames(ego) <- c("Category","ID","Term","Genes","adj_pval")
#attention：注意格式！
#善用class和str!
#这个包对格式要求很高！
ego$Genes <- str_replace_all(ego$Genes,"/",",")


circ <- circle_dat(ego, deg)
head(circ)
GOBar(subset(circ, category == 'BP'))     
GOBar(circ, display = 'multiple')      
GOBubble(circ, labels = 3)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
GOCircle(circ)
GOCircle(circ, nsub = c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047'))
GOCircle(circ, nsub = 5)

head(EC$genes)
EC$process
#注意数据类型
chord <- chord_dat(data = circ,genes = deg)
chord <- chord_dat(data = circ,process = ego$Term)
head(chord)
chord <- chord_dat(data = circ, genes = deg, process = ego$Term)

GOChord(data = chord, 
        title = "chord_plot",
        space = 0.02,
        limit = c(3,5),
        gene.order = 'logFC',
        gene.space = 0.25, 
        gene.size = 10,
        lfc.col = c("firebrick3","white","royalblue3"),
        ribbon.col = brewer.pal(length(ego$Term),"set3"),
        process.label = 18)


GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))


GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)


pdf("go_barplot.pdf",width=30,height=6)
dev.off()
