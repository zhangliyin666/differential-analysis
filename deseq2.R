count <- read.table("I:/ZCXrnaseq/genecounts.txt",header = T)
library(DESeq2)
library(tidyverse)
#refer ensembl gene id as standard
countdata <- count[,1:10]
treat <- countdata[,c(1,6:10)]
ctrl <- countdata[,1:5]
#There are two ways to merge them into one data set
#test1: by cbind
countdata_1 <- cbind(ctrl,treat)
countdata_1 <- countdata_1[,c(1:5,7:11)]
row.names(countdata_1) <- countdata_1[,1]
countdata_1 <- countdata_1[,-1]
head(countdata_1)
#test2: by merge
countdata_2 <- merge(ctrl,treat,by="gene_id")
row.names(countdata_2) <- countdata_2[,1]
countdata_2 <- countdata_2[,-1]
head(countdata_2)
###Here I use ensembl id and it needs to be transfer for convenient use.
###A problem is that when using gene name, an error appears and says rownames couldn't be duplicated values.
###If want to see the solution, see in another script "fpkm_to_tpm".

#grouping
condition <- factor(c(rep("control",4),rep("treat",5)), levels = c("control","treat"),ordered = F)
condition
table(condition)
coldata_1 <- data.frame(row.names=colnames(countdata_1), condition)
coldata_2 <- data.frame(row.names=colnames(countdata_2), condition)
coldata_2
table(coldata_2)

#to build dds object
dds_1 <- DESeqDataSetFromMatrix(countdata_1, coldata_1, design= ~ condition)
dds_2 <- DESeqDataSetFromMatrix(countdata_2, coldata_2, design= ~ condition)
#to filter low abundance data
dds_1 <- dds_1[rowSums(counts(dds_1)) > 1, ] 
dds_2 <- dds_2[rowSums(counts(dds_2)) > 1, ] 
#to calculate nomalized data
#colData(dds) adds column "sizeFcator", which is used to correct sequencing depth and component of library
dds_1 <- estimateSizeFactors(dds_1) 
dds_2 <- estimateSizeFactors(dds_2) 
#normalization
normalized_counts_1 <- counts(dds_1,normalized=T) 
normalized_counts_2 <- counts(dds_2,normalized=T)

#to set the order of comparison
#the first one is control group
dds_2$condition <- factor(dds_2$condition,levels = c("control","treat"))
dds_1$condition <- factor(dds_1$condition,levels = c("control","treat"))
#or
dds_2$condition <- relevel(dds_2$condition,ref = "control")
dds_1$condition <- relevel(dds_1$condition,ref = "control")

#differential analysis
dds_1 <- DESeq(dds_1)
dds_1
dds_2 <- DESeq(dds_2)
dds_2

#to get the result
res_1= results(dds_1)
res_1 = res_1[order(res_1$pvalue),]

res_2= results(dds_2)
res_2 = res_2[order(res_2$pvalue),]

#or
res <- results(dds, 
               contrast=c("condition","control","treat"),
               alpha = 0.05,#alpha: the adjusted p-value cutoff. 
               lfcThreshold = 0)
sum(res_1$pvalue < 0.05, na.rm=TRUE) 
sum(res_2$pvalue < 0.05 & res_2$log2FoldChange>1, na.rm=TRUE)
head(res_1)
head(res_2)

summary(res_1)
summary(res_2)
write.csv(res_1,file="All_results_1.csv")
write.csv(res_2,file="all_treat_vs_control.csv")

diff_gene_deseq2 <-subset(res_2, pvalue <= 0.05&abs(log2FoldChange) >2 )
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file = "diff_gene_deseq2.csv")
#annotation
library('biomaRt')
library("curl")
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff_gene_deseq2)
mmu_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)

head(mmu_symbols)
head(diff_gene_deseq2)
#get the annotated DEG
ensembl_gene_id<-rownames(diff_gene_deseq2)
diff_gene_deseq2<-cbind(ensembl_gene_id,diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
diff_name<-merge(diff_gene_deseq2,mmu_symbols,by="ensembl_gene_id")

write.csv(diff_name,file = "diff_name_deg.csv")
head(diff_name)

Ighg1 <- diff_name[diff_name$external_gene_name=="Ighg1",]
#boxplot
library(magrittr)
plotCounts(dds, gene="ENSMUSG00000076614", intgroup="condition", returnData=TRUE) %>% 
  ggplot(aes(condition, count)) + geom_boxplot(aes(fill=condition)) + scale_y_log10() + ggtitle("Akap8")
#point plot
d <- plotCounts(dds, gene="ENSMUSG00000076614", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(aes(color= condition),size= 4, position=position_jitter(w=0.5,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+ ggtitle("Ighg1")
#PCA plot by vst
vsdata <- vst(dds_2, blind=FALSE)
pcadata <- plotPCA(vsdata, intgroup="condition",returnData = T)
pcadata <- pcadata[order(pcadata$condition,decreasing=F),]
table(pcadata$condition)
#PCA after log2 tranferation
#to get expression set firstly
rld <- rlogTransformation(dds_2)
exprSet_new=assay(rld)
exprSet_new <- t(exprSet_new)
exprSet_new <- as.data.frame(exprSet_new)
exprSet_new <- cbind(exprSet_new,condition)
library(FactoMineR)
library(factoextra)
dat.pca <- PCA(exprSet_new[,-ncol(exprSet_new)], graph = F)
fviz_screeplot(dat.pca, addlabels = TRUE,
               ylim = c(0, 75))
plot(dat.pca,choix="ind")
var <- get_pca_var(dat.pca)
head(var$coord)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = exprSet_new$condition, # color by groups
             palette = c("#00AFBB", "#E7B800", '#CC00FF', '#FF0099'),
             addEllipses = TRUE, # Concentration ellipses # 是否圈起来
             legend.title = "Groups")

#count matrix heatmap
library("pheatmap")
select<-order(rowMeans(counts(dds_2, normalized = TRUE)),
              decreasing = TRUE)[1:50]#first 50 genes
df <- as.data.frame(colData(dds_2)[,c("condition","sizeFactor")])
ntd <- normTransform(dds_2)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=FALSE, annotation_col=df)
#sample-to-sample heatmap
sampleDists <- dist(t(assay(vsdata)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata@colData@rownames, vsdata$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



#enrichment by ClusterProfiler
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
sig.gene <- read.csv(file = "diff_name.csv")
#Attention the file format
#Make sure if the first column has colname
row.names(sig.gene) <- sig.gene$ensembl_gene_id
sig.gene <- sig.gene[,-1]
head(sig.gene)
#sometimes after input, some values turned into NA
#to find NAs and return the coordinates
which(is.na(sig.gene),arr.ind = T)
deg[is.na(sig.gene)] <- 0
#to set the up- or down- regulated genes
p_thred <- 0.05
logFC_thred <- 2
sig.gene$change=ifelse(sig.gene$pvalue>p_thred,'stable',
                       ifelse( sig.gene$log2FoldChange > logFC_thred,'UP',
                               ifelse( sig.gene$log2FoldChange < -logFC_thred,'DOWN','stable') )
)
table(sig.gene$change)
write.csv(sig.gene,file = "sig_gene.csv")
#assign the row names of sig.gene to gene (vector)
gene<-row.names(sig.gene)
head(gene)


#id transferation
gene.df<-bitr(gene, fromType = "ENSEMBL",
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Mm.eg.db)
head(gene.df)
gene_up <- sig.gene[sig.gene$change == 'UP','external_gene_name'] 
gene_down <- sig.gene[sig.gene$change == 'DOWN','external_gene_name'] 
gene_diff <- c(gene_up,gene_down)
#gene_all <- as.character(sig.gene[ ,'external_gene_name'] )


#GO enrichment
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
ego_cc@result$qvalue = 0
ego_cc@result$p.adjust = 0
dotplot(ego_cc,
        showCategory = 30,
        title="The GO_CC enrichment analysis of all DEGs ",
        color = 'pvalue')
ego_bp@result$qvalue = 0
ego_bp@result$p.adjust = 0
dotplot(ego_bp,
        showCategory = 30,
        title="The GO_BP enrichment analysis of all DEGs ",
        color = 'pvalue')
ego_mf@result$qvalue = 0
ego_mf@result$p.adjust = 0
dotplot(ego_mf,
        showCategory = 30,
        title="The GO_MF enrichment analysis of all DEGs ",
        color = 'pvalue')
meta <- ego_bp@result[grep(pattern="metaboli",ego_bp@result[,2]),]
#KEGG enrichment
kk<-enrichKEGG(gene      =gene.df$ENTREZID,
               organism = 'mmu',
               pvalueCutoff = 0.9)#0.05
#to convert the entrzid to symbol
kk=DOSE::setReadable(kk, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
kk[1:12]
kk@result$qvalue = 0
kk@result$p.adjust = 0
barplot(kk,showCategory = 30, title="The KEGG enrichment analysis of all DEGs",color = 'pvalue')
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs",color = 'pvalue')
#cnetplot(kk , categorySize="pvalue", foldChange=gene,colorEdge = TRUE)
cnetplot(kk,foldChange=gene, circular = TRUE, colorEdge = TRUE)
emapplot(kk,foldChange=gene, circular = TRUE, colorEdge = TRUE,color = "pvalue")
##GSEA 
sig.gene <- as.data.frame(diff_name)
gene_1=data.frame(sig.gene$ensembl_gene_id,sig.gene$log2FoldChange,stringsAsFactors=FALSE)
geneID=select(org.Mm.eg.db,keys=sig.gene$ensembl_gene_id,columns="ENTREZID",keytype="ENSEMBL")
geneID=na.omit(geneID)
colnames(gene_1) <- c("ENSEMBL", "log2FoldChange")
tmp=left_join(geneID,gene_1,by="ENSEMBL")
genelist=tmp$log2FoldChange
names(genelist)=tmp$ENTREZID
head(genelist)
genelist_sort=sort(genelist,decreasing = T)
head(genelist_sort)
#GO annotation
gsemf <- gseGO(genelist_sort,
               OrgDb = org.Mm.eg.db,
               keyType = "ENTREZID",
               ont="MF",
               pvalueCutoff = 0.9)
head(gsemf)
library(enrichplot)
gseaplot(gsemf, geneSetID=1,title = gsemf@result$Description[1])
gseaplot2(gsemf, geneSetID=1,title = gsemf@result$Description[1])
#KEGG annotation
kegg <- gseKEGG(genelist_sort ,
                organism     = 'mmu',
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.9,
                verbose      = FALSE)
head(kegg)

#pathway analysis
library("pathview")
library("gage")
library("gageData")
library("dplyr")
data("kegg.sets.mm")
data("sigmet.idx.mm")
kegg.sets.mm =  kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm,3)
#sig.gene<-read.csv(file="diff_gene_deseq2.csv")
#gene<-sig.gene[,1]
#head(gene)
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Mm.eg.db)
foldchanges = sig.gene$log2FoldChange
names(foldchanges)= gene.df$ENTREZID
head(foldchanges)
keggres = gage(foldchanges, gsets = kegg.sets.mm, same.dir = TRUE)
lapply(keggres, head)
keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
pics = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
#volcano plot 
library(ggplot2)
data <- as.data.frame(diff_name@listData)
data <- data[,c(1,3,6,8)]
head(data)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange > 1,'Up','Down'),
                     'Stable')
#method1
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = data$external_gene_name)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.000001),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential metabolites")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
p

library(ggrepel)
data$label=ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1,data$external_gene_name,"")
p+geom_text_repel(data = data, aes(x = data$log2FoldChange, 
                                   y = -log10(data$pvalue), 
                                   label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
#method2 
p <- ggplot(data = data, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue))) +
  geom_point(alpha=0.4, size=1, 
             aes(color=change)) +
  scale_color_manual(values=c("blue", "red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
#highlight the intersted genes
for_label <- data %>% 
  filter(abs(log2FoldChange) >3 & -log10(pvalue)> -log10(0.000001))
p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = external_gene_name),
    data = for_label,
    color="black"
  )


#Build the DEG heatmap using expression set
rld <- rlogTransformation(dds_2)
exprSet=assay(rld)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id_1<-row.names(exprSet)
mmu_symbols_1<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                      filters = 'ensembl_gene_id', values = my_ensembl_gene_id_1, mart = mart)
ensembl_gene_id_1<-rownames(exprSet)
exprSet<-cbind(ensembl_gene_id_1,exprSet)
colnames(exprSet)[1]<-c("ensembl_gene_id")
exprSet <- as.data.frame(exprSet)
exprSet_new<-merge(exprSet,mmu_symbols_1,by="ensembl_gene_id")
head(exprSet_new)#28367
write.csv(exprSet_new,file = "exprset_description.csv")
exprSet_new <- exprSet_new[,2:11]
b <- as.data.frame(exprSet_new$external_gene_name)
colnames(b) <- "symbol"
exprSet_new <- cbind(b,exprSet_new)
exprSet_new <- exprSet_new[,c(1:10)]
a <- exprSet_new
str(exprSet_new)
#exprSet_new <- a

#to get mean data because some observations appear two or more times
#ATTENTION: the type of data must be numeric
#You could write it into csv and then read it
write.csv(exprSet_new,'exprset.csv')
exprSet_new<- read.csv('exprset.csv',T)
exprSet_new <- exprSet_new[,-1]
#to take the average
data3 <- aggregate( . ~ symbol,data=exprSet_new, mean)
write.csv(data3,file = "exprset_filt.csv")
#or
data4 <- aggregate(exprSet_new[,2:10],list(exprSet_new$symbol),mean)

#DEG heatmap
FC <- sig.gene$log2FoldChange
names(FC) <- sig.gene$external_gene_name
class(FC)
FC
DEG_110 <- c(names(head(sort(FC),55)),names(tail(sort(FC),55)))
head(DEG_110)
rownames(data3) <- data3$symbol
data3 <- data3[,-1]
dat <- t(scale(t(data3[DEG_110,])))
dat[1:4,1:4]
group <- data.frame(group=condition)
rownames(group)=colnames(dat)
pheatmap(dat,show_colnames =T,show_rownames = T, 
         cluster_cols = F,cluster_rows = F,
         annotation_col=group,
         border_color = NA,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 

samples <- rep(c('control', 'treat'), c(4, 5))
heat <- Heatmap(dat, 
                col = colorRampPalette(c('navy', 'white', 'firebrick3'))(100), #to set the color from high to low expression
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #to set the height of legends
                show_row_names = T,  #not to show gene names
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('control' = '#00DAE0', 'treat' = '#FF9289')),  #to set the sample groups color
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 6),
                cluster_rows = F,
                cluster_columns = F)
heat
#show <- read.delim('show_name.txt', header = FALSE, check.names = FALSE)
#to show the interested genes
show <- as.data.frame(rownames(dat))
colnames(show)<-"V1"
show <- show[1:11,]
#!!!ATTENTION!!!
# "show" must be atomic vector but not recursive vector！
# "$" is for recursive and "[]" is for atomic!!
# is.recursive(show) and is.atomic(show) may help you
heat + rowAnnotation(link = anno_mark(at = which(rownames(dat) %in% show), 
                                      labels = show, labels_gp = gpar(fontsize = 5)))

