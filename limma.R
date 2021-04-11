library(dplyr)
library(edgeR)
library(limma)
library(edgeR)
count <- read.csv("H:/nash-analysis/gene_count.csv")
count <- filter(count,!duplicated(count$gene_name))
count1 <- count[,1:10]
treat1 <- count1[,c(1,8,2,3,7)]
ctrl1 <- count1[,c(1,4,10,6,9,5)]
count1 <- merge(ctrl1,treat1,by="gene_id")
#pre-processing
####Here I use ensembl id and it needs to be transfer for convenient use.
###A problem is that when using gene name, an error appears and says rownames couldn't be duplicated values.
###If want to see the solution, see in another script "fpkm_to_tpm".
#Here I just remove all duplicated data
#count1 <- round(count1)#round numbers
row.names(count1) <- count1[,1]
count1 <- count1[,-1]
head(count1)

#过滤低表达
d0 <- DGEList(count1)
d0 <- calcNormFactors(d0)
d0

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d0 <- d0[-drop,] 
#Groupinga nd Declare the order of contrast
condition <- factor(c(rep("control",5),rep("treat",4)), levels = c("control","treat"),ordered = F)
condition
table(condition)

design <- model.matrix(~0+condition)
colnames(design) = levels(factor(condition))
rownames(design) = colnames(condition)
design

#voom normalization
norm <- voom(d0, design, plot = TRUE)
#If the std is obviously high or if there is a rise at 0, it means that data includes many low counts
par(mfrow = c(2, 2))
boxplot(d0)
plotDensities(d0)
boxplot(norm$E)
plotDensities(norm$E)

#linear fitting
fit <- lmFit(norm, design, method = 'ls')

#To set the groups that ready to make contrast
#"1" means upregulation, "-1" means downregulation
contrast <- makeContrasts('treat-control', levels = design)
contrast

#The empirical Bayesian model was used to fit the standard deviation
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

qqt(fit2$t, df = fit2$df.prior+fit2$df.residual, pch = 16, cex = 0.2)
abline(0,1)

#The correction of p value and extraction of DEG
diff_gene <- topTable(fit2, number = Inf, adjust.method = 'fdr')
diff_gene <- subset(diff_gene,diff_gene$P.Value<=0.05)
head(diff_gene, 10)
write.csv(diff_gene, file = 'diff_gene_2.csv')

library("pheatmap")
exprset <- norm$E

#Relationship analysis
table(condition)
ac=data.frame(groups=condition)
rownames(ac)=colnames(exprset) #to set the column names of exprset as the row names of ac
pheatmap(cor(exprset),annotation_col = ac)

exprset=t(exprset) 
exprset=as.data.frame(exprset)
exprset=cbind(exprset,condition) #to add the group information to exprset

library("FactoMineR")
library("factoextra") 
exprset.pca <- PCA(exprset[,-ncol(exprset)], graph = FALSE)#The last column of exprset is group_list and pca analysis call a complete numeric matrix. So you should remove the last column.
(fviz_pca_ind(exprset.pca,geom.ind = "point",
              col.ind = exprset$condition,
              palette = c("#00AFBB", "#E7B800"),
              addEllipses = TRUE,
              legend.title = "Groups"
))
ggsave('all_samples_PCA.png')

#Check the exprset because we transpositted it
exprset <- norm$E
cg=names(tail(sort(apply(exprset,1,sd)),100))
pheatmap(exprset[cg,],show_colnames =F,show_rownames = F)
#Normalization
n=t(scale(t(exprset[cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
#add grouping information
ac=data.frame(g=condition)
rownames(ac)=colnames(n) 
(heatmap <- pheatmap(n,show_colnames =T,show_rownames = T,
                     annotation_col=ac,
                     cluster_rows = F,
                     cluster_cols = F))

deg <- diff_gene
head(deg)
save(deg,file = 'deg.Rdata')
#Volcano
df <- deg
colnames(deg)
df$v <- -log10(df$P.Value)
p_thred <- 0.05
logFC_thred <- 1
df$groups = ifelse(df$P.Value > p_thred, "stable", 
                   ifelse(df$logFC > logFC_thred, "up", 
                          ifelse(df$logFC < -logFC_thred, "down", 
                                 "stable")))
table(df$groups)

library(ggplot2)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_thred,3),
                    '\nThe number of up gene is ',nrow(df[df$groups =='up',]) ,
                    '\nThe number of down gene is ',nrow(df[df$groups =='down',])
)
p <- ggplot(data = df, aes(x = logFC, y = v)) +
  geom_point(alpha=0.4, size=1.75, 
             aes(color=groups)) +
  scale_color_manual(values=c("blue", "grey","red")) +
  geom_vline(xintercept=c(-logFC_thred,logFC_thred),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(p_thred),lty=4,col="black",lwd=0.8) +
  labs(x="logFC",y="-log10(P.value)")+
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  theme_bw()
p
#to highlight the interested genes
df$symbol <- row.names(df)
for_label <- df %>% 
  filter(abs(logFC) >0 & -log10(P.Value)> -log10(0.001))
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
#go enrichment
table(df$groups)
head(df)
df$SYMBOL=rownames(deg)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
gene.df <- bitr(unique(df$SYMBOL), fromType = "SYMBOL",
                toType = c( "ENTREZID","ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)
DEG=df
DEG=merge(DEG,gene.df,by='SYMBOL')
head(DEG)
table(DEG$groups)

save(DEG,file = 'anno_DEG.Rdata')
write.csv(DEG,file = "anno_DEG.csv")
######################################################

gene_up <- DEG[DEG$groups == 'up','ENTREZID'] 
gene_down <- DEG[DEG$groups == 'down','ENTREZID'] 
gene_diff <- c(gene_up,gene_down)
gene_all <- as.character(DEG[ ,'ENTREZID'] )

geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList,decreasing = T)


ego_cc_up<-enrichGO(gene       = gene_up,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENTREZID',
                    ont        = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
ego_cc_down<-enrichGO(gene       = gene_down,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
ego_cc_diff<-enrichGO(gene       = gene_diff,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)

ego_bp_up<-enrichGO(gene       = gene_up,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENTREZID',
                    ont        = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
ego_bp_down<-enrichGO(gene       = gene_down,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
ego_bp_diff<-enrichGO(gene       = gene_diff,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
ego_mf_up<-enrichGO(gene       = gene_up,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENTREZID',
                    ont        = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
ego_mf_down<-enrichGO(gene       = gene_down,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
ego_mf_diff<-enrichGO(gene       = gene_diff,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENTREZID',
                      ont        = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
dotplot(ego_cc_diff,
        showCategory = 20,
        title="The GO_CC_down enrichment analysis of all DEGs ")
#color = 'pvalue')
dotplot(ego_bp_diff,
        showCategory = 20,
        title="The GO_BP_diff enrichment analysis of all DEGs ")
#color = 'pvalue')
dotplot(ego_mf_diff,
        showCategory = 20,
        title="The GO_MF_diff enrichment analysis of all DEGs ")
#color = 'pvalue')

{
  
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Mm.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results.Rdata')
    
  }
  
  
  load(file = 'go_enrich_results.Rdata')
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}

#KEGG enrichment

kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'mmu',
                    universe     = gene_all,
                    pvalueCutoff = 0.05)
kk.down <- enrichKEGG(gene         = gene_down,
                      organism     = 'mmu',
                      universe     = gene_all,
                      pvalueCutoff = 0.05)
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'mmu',
                      universe     = gene_all,
                  pvalueCutoff = 0.05)

kk@result$qvalue = 0
kk@result$p.adjust = 0
kk.up=DOSE::setReadable(kk, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
barplot(kk.up,color = "pvalue")
cnetplot(kk.up, foldChange=gene_up,showCategory = 20, circular = TRUE, colorEdge = TRUE)
emapplot(kk,foldChange=gene, circular = TRUE, colorEdge = TRUE)

library("pathview")
library("gage")
library("gageData")
library("dplyr")
data("kegg.sets.mm")
data("sigmet.idx.mm")
kegg.sets.mm =  kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm,3)
foldchanges = DEG$logFC
names(foldchanges)= DEG$ENTREZID
head(foldchanges)
keggres = gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres, head)
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
