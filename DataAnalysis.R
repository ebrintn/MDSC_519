library(GEOquery)
library(Seurat)
library(rhdf5)
library(stringr)
library(Matrix)


#Read transcript data

#Control Data
control_data <- Read10X(data.dir = "SARS-CoV-2Dat/")
seurat_control <- CreateSeuratObject(counts = control_data, "control", min.cells = 10, min.features = 10)


#Type 2 diabetes data
T2D_data <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", sep = "\t")
T2D_dup <- T2D_data
T2D_data <- T2D_data[1:26179,]
n_occur <- data.frame(table(T2D_data[,1]))
duplicate_genes <- T2D_data[T2D_data[,1] %in% n_occur$Var1[n_occur$Freq > 1],]
duplicate_gene_groups <- split(rownames(duplicate_genes),duplicate_genes[,1])

for(group in duplicate_gene_groups){
  toRemove <- unlist(group)
  T2D_data[toRemove [1],] <- c(T2D_data[rownames(T2D_data) == toRemove [1],1:2],colSums(T2D_data[ rownames(T2D_data) %in% toRemove,3:7030]))
  T2D_data <- T2D_data[!(rownames(T2D_data) %in% toRemove[2:length(toRemove)]),]
  print(nrow(T2D_data))
}
rownames(T2D_data) <- T2D_data[,1]
cell_names_t2d <- read.delim("~/2020-21/MDSC519/Project/cell_names_t2d.txt", header=FALSE)
colnames(T2D_data) <- cell_names_t2d
counts_T2D_dat <- T2D_data[,3517:7030]
colnames(counts_T2D_dat) <- cell_names_t2d[2:3515]

write.csv(counts_T2D_dat, file = "refined_T2D_dat.csv")

counts_T2D_dat <- read.csv("refined_T2D_dat.csv", header = T, row.names = "X")

T2D_seurat <- CreateSeuratObject(counts = counts_T2D_dat, "T2D_dat")


#Add health data to the T2D data
cell_names_t2d <- read.delim("~/2020-21/MDSC519/Project/cell_names_t2d.txt", header=FALSE)
T2D_cell <- str_detect(cell_names_t2d[2:3515], "T2D")
health <- vector(mode="character", length=3514)
health[!T2D_cell] <- "control"
health[T2D_cell] <- "T2D"
health <- data.frame(health)
rownames(health) <- colnames(counts_T2D_dat)
T2D_seurat <- AddMetaData(T2D_seurat, health)



#Both T1D and T2D dataset
diabetes_data <- read.table("transcript_counts_final.tsv", sep = "\t", header = T, row.names = "target_id")
diabetes_data_refined <- diabetes_data[rowSums(diabetes_data)!=0,]


#Convert gene ids to gene names for diabetes data
gene_ids <- rownames(diabetes_data_refined)
gene_ids_conversion <- read.table("transcripts_to_genes.txt", row.names = "V1")
gene_ids_conversion <- gene_ids_conversion[rownames(gene_ids_conversion) %in% gene_ids,]
gene_ids_conversion <- gene_ids_conversion[gene_ids,]
rownames(gene_ids_conversion) <- gene_ids
gene_ids_conversion[is.na(gene_ids_conversion[,2]),] <- rownames(gene_ids_conversion[is.na(gene_ids_conversion[,2]),])
n_occur <- data.frame(table(gene_ids_conversion[,2]))
duplicate_genes <- gene_ids_conversion[gene_ids_conversion[,2] %in% n_occur$Var1[n_occur$Freq > 1] &! startsWith(gene_ids_conversion[,2], "ENST"),]
duplicate_gene_groups <- split(rownames(duplicate_genes),duplicate_genes[,2])

for(group in duplicate_gene_groups){
  toRemove <- unlist(group)
  diabetes_data_refined[rownames(diabetes_data_refined) == toRemove [1],] <- colSums(diabetes_data_refined[rownames(diabetes_data_refined) %in% toRemove,])
  diabetes_data_refined <- diabetes_data_refined[!(rownames(diabetes_data_refined)%in%toRemove[2:length(toRemove)]),]
  print(nrow(diabetes_data_refined))
}

write.csv(diabetes_data_refined, file = "refined_diabetes_dat.csv")

diabetes_data_refined <- read.csv("refined_diabetes_dat.csv", header = T, row.names = "X")

gene_ids <- rownames(diabetes_data_refined)
gene_ids_conversion <- read.table("transcripts_to_genes.txt", row.names = "V1")
gene_ids_conversion <- gene_ids_conversion[rownames(gene_ids_conversion) %in% gene_ids,]
gene_ids_conversion <- gene_ids_conversion[gene_ids,]
rownames(gene_ids_conversion) <- gene_ids
gene_ids_conversion[is.na(gene_ids_conversion$V2),] <- rownames(gene_ids_conversion[is.na(gene_ids_conversion$V2),])


diabetes_data_genenames <- diabetes_data_refined
rownames(diabetes_data_genenames) <- gene_ids_conversion[,2]

seurat_diabetes_dat <- CreateSeuratObject(counts = diabetes_data_genenames, "Diabetes Data")


#Add metadata to diabetes data
diabetes_metadata <- read.table("medata.txt", header = T, sep = ",", row.names = "Run")
diabetes_metadata <- diabetes_metadata[rownames(diabetes_metadata) %in% colnames(seurat_diabetes_dat),]
diabetes_metadata <- diabetes_metadata[, colnames(diabetes_metadata)]
seurat_diabetes_dat <- AddMetaData(seurat_diabetes_dat,diabetes_metadata)

toRemove <- colnames(seurat_diabetes_dat[,seurat_diabetes_dat$curated.cell.type == "dropped"])
diabetes_dat_filtered <- seurat_diabetes_dat[, !colnames(seurat_diabetes_dat) %in% toRemove]



#Perform QC and remove minimal transcripts

VlnPlot(diabetes_dat_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
diabetes_dat_filtered <- subset(diabetes_dat_filtered, subset = nFeature_RNA > 5000 & nFeature_RNA < 12500 & nCount_RNA < 3e6)

diabetes_dat_filtered <- NormalizeData(diabetes_dat_filtered, normalization.method = "LogNormalize", scale.factor = 1000)

VlnPlot(seurat_control, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seurat_control <- subset(seurat_control, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA < 3e4)

seurat_control <- NormalizeData(seurat_control, normalization.method = "LogNormalize", scale.factor = 10000)

VlnPlot(T2D_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
T2D_seurat <- subset(T2D_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA < 1e6)

T2D_seurat <- NormalizeData(T2D_seurat, normalization.method = "LogNormalize", scale.factor = 10000)




#Find variable features
diabetes_dat_filtered <- FindVariableFeatures(diabetes_dat_filtered)
top10 <- head(VariableFeatures(diabetes_dat_filtered), 10)
var_plot <- VariableFeaturePlot(diabetes_dat_filtered)
LabelPoints(plot = var_plot, points = top10, repel = TRUE)


seurat_control <- FindVariableFeatures(seurat_control)
top10 <- head(VariableFeatures(seurat_control), 10)
var_plot <- VariableFeaturePlot(seurat_control)
LabelPoints(plot = seurat_control, points = top10,repel = TRUE)


T2D_seurat <- FindVariableFeatures(T2D_seurat)
top10 <- head(VariableFeatures(T2D_seurat), 10)
var_plot <- VariableFeaturePlot(T2D_seurat)
LabelPoints(plot = T2D_seurat, points = top10,repel = TRUE)


#Integrate datasets
features <- SelectIntegrationFeatures(object.list = c(diabetes_dat_filtered, seurat_control, T2D_seurat))
diabetes.anchors <- FindIntegrationAnchors(object.list = c(diabetes_dat_filtered, seurat_control, T2D_seurat), anchor.features = features)
diabetes_combined <- IntegrateData(anchorset = diabetes.anchors)



#Scale Data in 3 datasets


all.genes <- rownames(diabetes_combined)
diabetes_combined<- ScaleData(diabetes_combined, features = all.genes)




#Run PCA and find number of principal components

diabetes_combined <- RunPCA(diabetes_combined)
DimPlot(diabetes_combined, reduction = "pca")
DimHeatmap(diabetes_combined, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(diabetes_combined, ndims = 50)



#Find clusters and nearest neighbours for 3 datasets

diabetes_combined <- FindNeighbors(diabetes_combined, dims = 1:25)
diabetes_combined <- FindClusters(diabetes_combined, resolution = 0.5)

diabetes_combined <- RunUMAP(diabetes_combined, dims = 1:25)



#Plot UMAP of combined datasets by cell type and health
diabetes_combined$health[is.na(diabetes_combined$health)] <- "control"

Idents(diabetes_combined) <- "seurat_clusters"
DimPlot(diabetes_combined, reduction = "umap")

Idents(diabetes_combined) <- "curated.cell.type"
DimPlot(diabetes_combined, reduction = "umap")


#Find Cell Identities

acinar <- c("PRSS1", "PNLIP", "CTRC", "PNLIPRP1")
alpha <- c("GCG", "FXYD5", "CHGA", "PCSK2")
beta <- c("SLC30A8", "HOPX", "MAFA", "INS")
delta <- c("FRZB", "MS4A8", "CASR", "BCHE")
ductal <- c("PDLIM3", "CLDN1", "LGALS4", "WFDC2")
epsilon <- c("SPTSSB", "APOH", "CORIN", "ADAMTS6")
pp <- c("SCGN", "SCGB2A1", "ZNF503", "PPY")
stellate <- c("TIMP1", "NDUFA4L2", "SPON2", "GEM")

FeaturePlot(diabetes_combined, features = acinar)
FeaturePlot(diabetes_combined, features = alpha)
FeaturePlot(diabetes_combined, features = beta)
FeaturePlot(diabetes_combined, features = delta)
FeaturePlot(diabetes_combined, features = ductal)
FeaturePlot(diabetes_combined, features = epsilon)
FeaturePlot(diabetes_combined, features = pp)
FeaturePlot(diabetes_combined, features = stellate)



DoHeatmap(diabetes_combined, features = acinar) 
DoHeatmap(diabetes_combined, features = alpha) 
DoHeatmap(diabetes_combined, features = beta)
DoHeatmap(diabetes_combined, features = delta)
DoHeatmap(diabetes_combined, features = ductal)
DoHeatmap(diabetes_combined, features = epsilon)
DoHeatmap(diabetes_combined, features = pp)
DoHeatmap(diabetes_combined, features = stellate)

DoHeatmap(diabetes_combined, features = c(acinar,alpha,beta,delta,ductal,epsilon,pp,stellate))


new_cluster_ids <- c("Acinar", "Ductal","Alpha","Beta","PP","PP","Acinar","PP","Stellate","Ductal","Delta","Alpha","Unknown","Alpha","Beta","Alpha","Stellate")
names(new_cluster_ids) <- levels(diabetes_combined)
diabetes_combined <- RenameIdents(diabetes_combined, new_cluster_ids)
curated.cell.type <- Idents(diabetes_combined)
T2D_seurat <- AddMetaData(diabetes_combined, Idents(diabetes_combined))
DimPlot(diabetes_combined, reduction = "umap", label = TRUE, pt.size = 0.5)  + NoLegend()


#Analyze ACE2 and TMPRSS2 expression
VlnPlot(diabetes_combined, "ACE2", split.by = "health") + aes(color = diabetes_combined$health) 
VlnPlot(diabetes_combined, "ACE2") + aes(color = Idents(diabetes_combined)) + NoLegend()

VlnPlot(diabetes_combined, "TMPRSS2" , split.by = "health") + aes(color = diabetes_combined$health) 
VlnPlot(diabetes_combined, "TMPRSS2") + aes(color = Idents(diabetes_combined)) + NoLegend()

FeaturePlot(diabetes_combined, features = c("ACE2"), min.cutoff = -0.001, max.cutoff = 0.001, split.by = "health", keep.scale = "all", cols = c("Red","Blue", "Green"), pt.size = 0.05)
FeaturePlot(diabetes_combined, features = c("TMPRSS2"), min.cutoff = 0, max.cutoff = 0.1, split.by = "health", keep.scale = "all", cols = c("Blue", "Green"), pt.size = 0.05)
DoHeatmap(diabetes_combined, features = c("ACE2")) 


#Look at only beta cells
VlnPlot(diabetes_combined, "ACE2", split.by = "health", idents = "Beta") + aes(color = diabetes_combined$health[Idents(diabetes_combined)=="Beta"])+NoLegend()
VlnPlot(diabetes_combined, "TMPRSS2", split.by = "health", idents = "Beta") + aes(color = diabetes_combined$health[Idents(diabetes_combined)=="Beta"])+NoLegend()


#Look at only alpha cells
VlnPlot(diabetes_combined, "ACE2", split.by = "health", idents = "Alpha") + aes(color = diabetes_combined$health[Idents(diabetes_combined)=="Alpha"])+NoLegend()
VlnPlot(diabetes_combined, "TMPRSS2", split.by = "health", idents = "Alpha") + aes(color = diabetes_combined$health[Idents(diabetes_combined)=="Alpha"])+NoLegend()
