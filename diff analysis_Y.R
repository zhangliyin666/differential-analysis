data = read.csv("gene_countshm3.csv") 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
data <- read.csv('gene_countshm3.csv')
countdata <- data[,1:7]
#gene <- data[,8]
#countdata <- cbind(gene,countdata)
row.names(countdata) <- countdata[,1]
countdata <- countdata[,-1]
head(countdata)
condition <- factor(c(rep("control",3),rep("kd",3)), levels = c("control","kd"))
condition
coldata <- data.frame(row.names=colnames(countdata), condition)
coldata

dds <- DESeqDataSetFromMatrix(countdata, coldata, design= ~ condition)
dds <- DESeq(dds)
dds

res = results(dds, contrast=c("condition", "control", "kd"))
head(res)
summary(res)
res$padj[is.na(res$padj)] <- 1
write.csv(res,file="All_results.csv")

table(res$pvalue<0.05)
table(res$log2FoldChange<-1)
table(res$log2FoldChange>1)
fivenum(res$log2FoldChange<0)
diff_gene_deseq2$pvalue
table(res$padj<0.05 & abs(res$log2FoldChange>1))
diff_gene_deseq2 <-subset(res,((padj < 0.05)) & (log2FoldChange > 1 | log2FoldChange < -1))
#diff_gene_deseq4 <-subset(diff_gene_deseq3,(log2FoldChange > 0.5 | log2FoldChange < -0.5))
#diff_gene_deseq3 <-subset(diff_gene_deseq2,((padj < 0.05)|(padj == NA)))
#head(diff_gene_deseq3)
dim(diff_gene_deseq2)
summary(diff_gene_deseq2)
log2(1.2)
table(diff_gene_deseq2$log2FoldChange<-1)
table(diff_gene_deseq2$log2FoldChange>1)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_treat_vs_control.csv")

library('biomaRt')
library("curl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff_gene_deseq2)
hsa_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(hsa_symbols)

head(diff_gene_deseq2)
head(hsa_symbols)

ensembl_gene_id<-rownames(diff_gene_deseq2)
diff_gene_deseq2<-cbind(ensembl_gene_id,diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
diff_name<-merge(diff_gene_deseq2,hsa_symbols,by="ensembl_gene_id")
head(diff_name)
METTL3 <- diff_name[diff_name$external_gene_name=="METTL3",]
METTL3



#BiocManager::install("apeglm")
#res.shrink <- lfcShrink(dds, contrast = c("condition","kd","control"), res=res,type = 'normal')
#plotMA(res.shrink, ylim = c(-5,5))
#topGene <- rownames(res)[which.min(res$padj)]
#with(res[topGene, ], {
 # points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  #text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
#})

#plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
#plotCounts(dds, gene="ENSG00000156508", intgroup="condition", returnData=FALSE)

#plotCounts(dds, gene="ENSG00000000971", intgroup="condition", returnData=TRUE) %>% 
 # ggplot(aes(condition, count)) + geom_boxplot(aes(fill=condition)) + scale_y_log10() + ggtitle("CFH")

#vsdata <- vst(dds, blind=FALSE)
#plotPCA(vsdata, intgroup="condition")

library("pheatmap")
select<-order(rowMeans(counts(dds, normalized = TRUE)),
              decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsdata)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$condition, vsdata$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
sig.gene<-read.csv(file="DEG_treat_vs_control.csv")
head(sig.gene)
gene<-sig.gene[,1]
head(gene)
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)

ego_cc<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
               
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05
                 )
ego_bp<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf@result$qvalue = 0
ego_mf@result$p.adjust = 0
clusterProfiler::dotplot(ego_mf,showCategory = 20,title="The GO_MF enrichment analysis of all DEGs ",color = 'pvalue')
dotplot(ego_bp,showCategory = 20,title="The GO_BP enrichment analysis of all DEGs ")
dotplot(ego_mf,showCategory = 20,title="The GO_MF enrichment analysis of all DEGs ")
?dotplot
library(stringr)
kk<-enrichKEGG(gene      =gene.df$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.05)#0.05
kk[1:30]
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")
dotplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")
#cnetplot(kk , categorySize="pvalue", foldChange=gene,colorEdge = TRUE)
cnetplot(kk, foldChange=gene, circular = TRUE, colorEdge = TRUE)
emapplot(kk)
class(tmp)
str(tmp)

tmp=read.csv('All_results.csv')
tmp = tmp[abs(tmp$log2FoldChange)>0.5,]
tmp = tmp[tmp$pvalue<0.05,]
tmp = tmp[!duplicated(tmp$X),]
upgenes = tmp[tmp$log2FoldChange>0,]$X
downgenes = tmp[tmp$log2FoldChange<0,]$X



tmp2 = tmp[1:1749,]
res@listData[1:3,1:3]
tmp2=tmp2[abs(tmp2$log2FoldChange)>0.5,]
tmp2=tmp2[tmp2$pvalue<0.05,]
tmp2= as.data.frame.matrix(tmp)
tmp2=tmp2[!duplicated(tmp2$X),]
upgenes = tmp2[tmp2$log2FoldChange>0,]$X
downgenes = tmp2[tmp2$log2FoldChange<0,]$X
upgenes%in%downgenes
table(ego_cc@result$ID)
fivenum(ego_cc@result$pvalue)
tmp2 = tmp2[-1*is.na(tmp2$baseMean),]
tail(tmp2)
tmp = subset(tmp,)
tmp= as.matrix.data.frame(tmp)
tmp = res[res$log2FoldChange!=NA,]
head(tmp)

genelist <- sig.gene$log2FoldChange
names(genelist) <- sig.gene[,1]
genelist <- sort(genelist, decreasing = TRUE)
gsemf <- gseGO(genelist,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont="BP",
               pvalueCutoff = 0.9)
head(gsemf)
gseaplot(gsemf, geneSetID="GO:0001935")

BiocManager::install("pathview")
BiocManager::install("gage")
BiocManager::install("gageData")
library(pathview)
library(gage)
library(gageData)
library(dplyr)
#library(clusterProfiler)
#library(DOSE)
#library(stringr)
#library(org.Mm.eg.db)

data("kegg.sets.hs")
data("sigmet.idx.hs")
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
sig.gene<-read.csv(file="DEG_treat_vs_control.csv")
gene.df<-bitr(gene, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(sig.gene)

foldchanges = sig.gene$log2FoldChange
names(foldchanges)= gene.df$ENTREZID
head(foldchanges)

keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
lapply(keggres, head)

keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()
keggrespathways

keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))


threshold<-as.factor(
  (sig.gene$log2FoldChange>1|sig.gene$log2FoldChange<(-1) & 
     sig.gene$pvalue<0.05 ))
library(ggplot2)
ggplot(sig.gene,aes(x= log2FoldChange,
               y= -1*log10(pvalue),colour=threshold))+xlab("log2 fold-change")+ylab("-log10 p-value")+geom_point() 
