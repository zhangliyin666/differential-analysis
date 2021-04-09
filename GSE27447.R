#step1-读入数据并获取分组信息
rm(list=ls())
options(stringsAsFactors = F)

library(AnnoProbe) # 生信技能树出品的GEO数据下载利器
library(Biobase)
gset <- geoChina("GSE27447")
gse27447 <- gset[[1]]
exprSet <- exprs(gse27447)
boxplot(exprSet,las=2)

#根据boxplot可以得出，原始数据需要log2处理。
exprSet <- log2(exprSet+1)
boxplot(exprSet,las=2)

#行名是探针，列名是样本，中间的数据是某样本中某探针的表达量。
exprSet[1:4,1:4]

# checkGPL()结果是TRUE说明AnnoProbe包中存在"GPL6244" 平台数据，于是可以使用另外两个超级厉害的函数idmap()和filterEM()，得到探针对应的基因名，然后把表达矩阵的探针名转换为基因名。
gse27447@annotation
checkGPL(gse27447@annotation)

#获取临床信息，从中进一步获取分组信息
ids <- idmap(gse27447@annotation)
dat <- filterEM(exprSet,ids)
dim(dat)
#[1] 18837    19
dat <- dat[order(rownames(dat)),]
pd <- pData(gse27447)#查看样本分组信息
library(stringr)
group_list=str_split(pd$title,' ',simplify = T)[,1]#按title分组
table(group_list)
#所以是5个三阴乳腺癌样本，14个非三阴乳腺癌样本。

#保存
save(dat,group_list,file = 'step1-output.Rdata')

#step2-检查表达矩阵
rm(list = ls())  
options(stringsAsFactors = F)
load('step1-output.Rdata')
#每一次都要检查
table(group_list)
dat[1:4,1:4]
#样本相关性分析
dim(dat)
ac=data.frame(groups=group_list)
rownames(ac)=colnames(dat) #把ac的行名给到n的列名，即对每一个探针标记上分组信息
pheatmap(cor(dat),annotation_col = ac)
#主成分分析
#对样本做主成分分析，要求行名是样本名，列名是探针（基因）名，所以需要对表达矩阵进行转置
dat=t(dat) #转置
dat=as.data.frame(dat)
dat=cbind(dat,group_list) ##给表达矩阵加上分组信息

library("FactoMineR")
library("factoextra") 
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)#dat最后一列是group_list。pca分析需要一个纯数值矩阵，所以将dat最后一列去掉以后赋值给dat.pca
(fviz_pca_ind(dat.pca,geom.ind = "point",
              col.ind = dat$group_list,
              palette = c("#00AFBB", "#E7B800"),
              addEllipses = TRUE,
              legend.title = "Groups"
))
ggsave('all_samples_PCA.png')
#取sd最大的1000个基因画热图
rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'step1-output.Rdata')
dat[1:4,1:4] 
##使用 apply()函数计算获取表达矩阵dat每行（即每个基因表达量）的方差，从小到大排序后，取最大的1000个方差，获取其对应的基因名，赋值给变量cg
##然后就可以对这1000个基因画热图
cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =T,show_rownames = T)
#对dat[cg,]归一化
n=t(scale(t(dat[cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
#给样本添加分组信息
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个样本标记上分组信息
(heatmap <- pheatmap(n,show_colnames =F,show_rownames = F,
                     annotation_col=ac))
library(ggplot2)
ggsave(filename = 'heatmap_top1000_sd.png',heatmap)

#step3-差异分析
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output.Rdata') ##载入数据
dat[1:4,1:4] 
group_list[group_list=='non-triple'] <- 'non_triple'
table(group_list) 
boxplot(unlist(dat[1,])~group_list)  #按照group_list分组画箱线图
bp=function(g){         #定义一个函数g，函数为{}里的内容
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
#随便检查几个基因的数据分布
p1 <- bp(as.numeric(dat[1,]))
p2 <- bp(as.numeric(dat[2,]))
p3 <- bp(as.numeric(dat[3,]))
p4 <- bp(as.numeric(dat[4,]))
library(patchwork)
(bp_4 <- p1|p2|p3|p4)
#差异分析正式开始
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
exprSet=dat
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts("triple-non_triple",levels = design)
contrast.matrix ##这个矩阵声明，我们要把 triple 组跟 non_triple 进行差异分析比较
#自定义函数DEG，功能是做差异分析。
DEG <- function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2) 
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  return(nrDEG)
}
#赋给deg
deg <- DEG(exprSet,design,contrast.matrix)
head(deg)
#保存
save(deg,file = 'deg.Rdata')
#手动检查deg前几个基因表达情况
bp(as.numeric(dat[rownames(dat)=='CR2',]))
#差异分析结果可视化
#volcano plot
#设置P.Value和logFC的阈值，将基因分为up，down，stable三类
df <- deg
colnames(deg)
df$v <- -log10(df$P.Value)
p_thred <- 0.05
logFC_thred <- 1.5
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
print(p)
#heatmap
load(file = 'step1-output.Rdata')
dat[1:4,1:4]
table(group_list)
x=deg$logFC #deg取logFC这列并将其重新赋值给x
names(x)=rownames(deg) #deg取probe_id这列，并将其作为名字给x
cg=c(names(head(sort(x),100)),#对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
     names(tail(sort(x),100)))

n=t(scale(t(dat[cg,])))
n[n>2]=2
n[n< -2]= -2

ac=data.frame(groups=group_list)
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac) 

#step4-GO/KEGG数据库注释
#整理数据
#设置P.Value和logFC的阈值，将基因分为UP，DOWN，stable3组，deg矩阵新建一列g储存基因分组信息。
#不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'deg.Rdata')
head(deg)
p_thred <- 0.05
logFC_thred <- 1.5
deg$g=ifelse(deg$P.Value>p_thred,'stable',
             ifelse( deg$logFC > logFC_thred,'UP',
                     ifelse( deg$logFC < -logFC_thred,'DOWN','stable') )
)
table(deg$g)
head(deg)
#富集分析需要基因的ENTREZID，使用Y叔神作clusterProfiler包里的bitr()函数获取对应关系，deg矩阵转存为DEG，并加上ENTREZID列
deg$SYMBOL=rownames(deg)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$SYMBOL), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)
#保存
save(DEG,file = 'anno_DEG.Rdata')
#分别取出上调基因和下调基因，合并为差异基因
gene_up <- DEG[DEG$g == 'UP','ENTREZID'] 
gene_down <- DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff <- c(gene_up,gene_down)
gene_all <- as.character(DEG[ ,'ENTREZID'] )

geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList,decreasing = T)

#KEGG
library(ggplot2)
library(clusterProfiler)
kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
head(kk.up)[,1:6]
#KEGG分析上调基因集结果可视化
#上调基因所属信号通路（气泡图）
dotplot(kk.up)
kk.up=DOSE::setReadable(kk.up, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk.up 

#查看第一个结果hsa04640的信号通路示意图
browseKEGG(kk.up, 'hsa04640')

#上调基因所属信号通路（条带图）
(gg_barplot <- barplot(kk.up,showCategory=20))

#通路与基因之间的关系可视化
cnetplot(kk.up , categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(kk.up, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

#通路与通路之间的关系图
emapplot(kk.up)

#通路与基因之间的关系(热图)
heatplot(kk.up)

#同样的方法计算下调基因集KEGG
kk.down <- enrichKEGG(gene         =  gene_down,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
head(kk.down)[,1:6]
dotplot(kk.down )
kk.down=DOSE::setReadable(kk.down, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk.down 
#查看第一个结果hsa04640的信号通路示意图
browseKEGG(kk.up, 'hsa04915')
#上调基因所属信号通路（条带图）
(gg_barplot <- barplot(kk.down,showCategory=20))

#通路与基因之间的关系可视化
cnetplot(kk.down , categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(kk.down, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

#通路与通路之间的关系图
emapplot(kk.down)

#通路与基因之间的关系(热图)
heatplot(kk.down)

# 差异基因集
kk.diff <- enrichKEGG(gene         = gene_diff,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
dotplot(kk.diff );ggsave('kk.diff.dotplot.png')

# Pathway Enrichment
kegg_diff_dt <- as.data.frame(kk.diff)
kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = 'kegg_up_down.png')
#这张图其实就是上调基因和下调基因pathway条带图合并，可以看到Oxytocin signaling pathway同时存在于上调基因和下调基因pathway，刚好跟差异基因集的分析结果一样。

#GSEA
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))

down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1

g_kegg=kegg_plot(up_kegg,down_kegg)
print(g_kegg)
ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')

#GO
g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

go_enrich_results <- lapply( g_list , function(gene) {
  lapply( c('BP','MF','CC') , function(ont) {
    cat(paste('Now process ',ont ))
    ego <- enrichGO(gene          = gene,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
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

load(file = 'go_enrich_results.Rdata')
n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1500,height = 1500)
    dotplot(go_enrich_results[[i]][[j]],title=paste0('dotplot_',n1[i],'_',n2[j])) %>% print()#记得勾选magrittr包
    dev.off()
  }
}
  

  
  

#step5-anno-GSEA
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

### 对 MigDB中的全部基因集 做GSEA分析。
# http://www.bio-info-trainee.com/2105.html
# http://www.bio-info-trainee.com/2102.html 
{
  load(file = 'anno_DEG.Rdata')
  geneList=DEG$logFC
  names(geneList)=DEG$symbol
  geneList=sort(geneList,decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  d='../MsigDB/symbols'
  gmts <- list.files(d,pattern = 'all')
  gmts
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
  ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
  f='gsea_results.Rdata'
  if(!file.exists(f)){
    gsea_results <- lapply(gmts, function(gmtfile){
      # gmtfile=gmts[2]
      geneset <- read.gmt(file.path(d,gmtfile)) 
      print(paste0('Now process the ',gmtfile))
      egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
      head(egmt)
      # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
      
      return(egmt)
    })
    # 上面的代码耗时，所以保存结果到本地文件
    save(gsea_results,file = f)
  }
  load(file = f)
  #提取gsea结果，熟悉这个对象
  gsea_results_list<- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  gseaplot(gsea_results[[2]],'KEGG_CELL_CYCLE') 
  gseaplot(gsea_results[[2]],'FARMER_BREAST_CANCER_CLUSTER_6') 
  gseaplot(gsea_results[[5]],'GO_CONDENSED_CHROMOSOME_OUTER_KINETOCHORE') 
  
}

#step6-anno-GSVA
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

### 对 MigDB中的全部基因集 做GSVA分析。
## 还有ssGSEA, PGSEA
{
  load(file = 'step1-output.Rdata')
  # 每次都要检测数据
  dat[1:4,1:4]  
  
  X=dat
  table(group_list)
  ## Molecular Signatures Database (MSigDb) 
  d='../MSigDB/symbols/'
  gmts=list.files(d,pattern = 'all')
  gmts
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSVA) # BiocManager::install('GSVA')
  
  if(T){
    es_max <- lapply(gmts, function(gmtfile){ 
      #gmtfile=gmts[8];gmtfile
      geneset <- getGmt(file.path(d,gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff=FALSE, verbose=FALSE, 
                     parallel.sz=1)
      return(es.max)
    })
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)
    es_deg <- lapply(es_max, function(es.max){
      table(group_list)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list))
      colnames(design)=levels(factor(group_list))
      rownames(design)=colnames(es.max)
      design
      library(limma)
      contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),
                                     levels = design)
      contrast.matrix<-makeContrasts("Tumor-Normal",
                                     levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
      
      deg = function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max,design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value=adjPvalueCutoff)
        summary(res)
        tempOutput = topTable(fit2, coef=1, n=Inf)
        nrDEG = na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re = deg(es.max,design,contrast.matrix)
      nrDEG=re
      head(nrDEG) 
      return(nrDEG)
    })
  } 
  
  gmts
  
  save(es_max,es_deg,file='gsva_msigdb.Rdata')
  
  load(file='gsva_msigdb.Rdata')
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=2
    print(i)
    dat=es_max[[i]]
    df=es_deg[[i]]
    df=df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n=rownames(df)
      dat=dat[match(n,rownames(dat)),]
      ac=data.frame(g=group_list)
      rownames(ac)=colnames(dat)
      rownames(dat)=substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'.pdf'))
      
    }
  })
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df=do.call(rbind ,es_deg)
  es_matrix=do.call(rbind ,es_max)
  df=df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df,file = 'GSVA_DEG.csv')
}
