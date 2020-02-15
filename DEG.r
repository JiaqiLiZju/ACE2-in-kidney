library(Seurat)
library(ggplot2)
library(ggrepel)

data <- read.table("dge.csv", sep=',', header=T, row.names="index")
kid <- CreateSeuratObject(counts =t(data)) 

PT_ACE2_label <- read.table("dge_obs.csv", sep=',', header=T, row.names="index")$PT_ACE
kid@meta.data$PT_ACE2_label <- as.factor(PT_ACE2_label)
# list options for groups to perform differential expression on
levels(kid)
# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
marker <- FindMarkers(kid, ident.1="PT_ACE2_pos", ident.2="PT_ACE2_neg", group.by='PT_ACE2_label', min.pct = 0.2)
hvg_gene <- marker
summary(hvg_gene)
# 
hvg_gene$threshold = factor(ifelse(hvg_gene$avg_logFC <= -1 | hvg_gene$avg_logFC >= 1, 
                                   ifelse(hvg_gene$avg_logFC >= 0,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
table(hvg_gene$threshold)

hvg_gene$gene<-rownames(hvg_gene)

hvg_gene<-cbind(hvg_gene,-log10(hvg_gene$p_val_adj))
colnames(hvg_gene)[ncol(hvg_gene)]<-"log_p_val_adj"
#hvg_gene <- hvg_gene[-1,]

p <- ggplot(hvg_gene,aes(x=avg_logFC,y=log_p_val_adj,color=threshold))+
  geom_point(size=1.5)+
  scale_color_manual(values=c("#FF0033","#0099FF","#000000"))+#确定点的颜色
  geom_text_repel(
    data = hvg_gene[hvg_gene$avg_logFC <= -1 | hvg_gene$avg_logFC >= 1,],
    aes(label = gene),
    size = 5,
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  theme_bw()+#修改图片背景
  theme(
    legend.title = element_blank()#不显示图例标题
  )+
  ylab('-log10 (p_val_adj)')+#修改y轴名称
  xlab('avg_logFC')+#修改x轴名称
  # geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  # geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+#添加竖线padj<0.05
  xlim(-3,3)+
  ylim(0,150)
 
pdf("gene_reg.pdf");p;dev.off()
write.csv(marker, file = "marker.csv")
save.image("workspace.rdata")