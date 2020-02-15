library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)

pdf("7.Dim.pdf", width=10, height=7);p;dev.off()

#Kidney data loading
K1.data <- Read10X(data.dir = "../data/GSE131685/GSM4145204_kidney1/")
K1 <- CreateSeuratObject(counts = K1.data, project = "kidney1", min.cells = 8, min.features = 200)
K2.data <- Read10X(data.dir = "../data/GSE131685/GSM4145205_kidney2/")
K2 <- CreateSeuratObject(counts = K2.data, project = "kidney2", min.cells = 6, min.features = 200)
K3.data <- Read10X(data.dir = "../data/GSE131685/GSM4145206_kidney3/")
K3 <- CreateSeuratObject(counts = K3.data, project = "kidney3", min.cells = 10, min.features = 200)
kid <- merge(x = K1, y = list(K2, K3))

#quality control
kid[["percent.mt"]] <- PercentageFeatureSet(kid, pattern = "^MT-")
VlnPlot(kid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
kid <- subset(kid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30)
kid <- NormalizeData(kid, normalization.method = "LogNormalize", scale.factor = 10000)
kid <- NormalizeData(kid)
kid <- FindVariableFeatures(kid, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kid), 10)
plot1 <- VariableFeaturePlot(kid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kid <- CellCycleScoring(kid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.genes <- rownames(kid)
kid <- ScaleData(kid, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

#Eliminate batch effects with harmony and cell classification
kid <- RunPCA(kid, pc.genes = kid@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
kid <- kid %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(kid, 'harmony')
harmony_embeddings[1:5, 1:5]

kid <- kid %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.25) %>% 
    identity()

new.cluster.ids <- c(1,2,3,4,5,6,7,8,9)
new.cluster.ids <- c("Proximal convoluted tubule",
    "Proximal tubule",
    "Proximal straight tubule",
    "NK cell and T cell",
    "Monocytes and B cells",
    "Glomerular parietal epithelial cells",
    "Distal tubule",
    "Collecting duct principal cells",
    "Collecting duct intercalated cells")
names(new.cluster.ids) <- levels(kid)
kid <- RenameIdents(kid, new.cluster.ids)

DimPlot(kid, reduction = "umap", label=TRUE, label.size = 4)
FeaturePlot(kid, features = c('ACE2', "CUBN", "LRP2"))
VlnPlot(kid, features = c('ACE2'))

# select cells
PT_select = kid@meta.data$seurat_clusters%in%c(0,1,2)
PT_ACE2_pos_sel = PT_select & kid@assays$RNA['ACE2']>0
PT_ACE2_neg_sel = PT_select & kid@assays$RNA['ACE2']<=0
l=length(kid@meta.data$seurat_clusters)
PT_ACE2_label = c(rep('NOT_PT', l))
PT_ACE2_label[PT_ACE2_pos_sel[1:l]] <- 'PT_ACE2_pos'
PT_ACE2_label[PT_ACE2_neg_sel[1:l]] <- 'PT_ACE2_neg'
table(PT_ACE2_label)

# updata ident
kid@meta.data$PT_ACE2_label <- as.factor(PT_ACE2_label)
kid@active.ident <- as.factor(PT_ACE2_label)
p <- DimPlot(kid, reduction = "umap", cols=c('green','grey','red'),label=FALSE, label.size = 4)

# restore ident
kid@active.ident <- as.factor(kid@meta.data$seurat_clusters)

#Calculating differentially expressed genes (DEGs) and Save rds file
kid.markers <- FindAllMarkers(kid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(kid.markers, sep=",", file="markers.csv")
saveRDS(kid, file="data.rds")
save.image("workspace.rdata")