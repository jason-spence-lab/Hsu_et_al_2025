library(Seurat)
library(ggplot2)
Dotplot_Zhiwei_Version <- function(seurat_object, gene_list, grouping) {
  DotPlot(seurat_object, features = gene_list,  dot.scale = 20, scale = TRUE, group.by = grouping) + RotatedAxis() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
    scale_colour_distiller(palette = "YlGnBu", direction = 1) +
    guides(size=guide_legend(title = "Percent Expressed",override.aes=list(shape=21, colour="black", fill="black"))) +
    labs(y=NULL, x= NULL) +
    guides(color = guide_colourbar(title = "Scaled Expression", ticks = TRUE, frame.colour = "black")) +
    theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())+
    theme(panel.border=element_rect(colour="black",fill=NA,size=0.8))
}
#load cellbender corrected objects, requires ReadCB_h5 function to load
library(Matrix)
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- "matrix"
  output <- list()
  
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # CellRanger version 3
    message('CellRanger version 3+ format H5')
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    message('CellRanger version 2 format H5')
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    
    if (unique.features) {
      features <- make.unique(names = features)
    }
    
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'CsparseMatrix')
    
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    
    output[[genome]] <- sparse.mat
  }
  
  infile$close_all()
  
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else {
    return(output)
  }
}
counts.trachea.one <- ReadCB_h5("~/scRNAseq/Data/Multiome/Ansley and Peggy Proximal to Distal scRNA-seq/6519-AC/Sample_6519-AC-1/Cellbender v03/output_filtered.h5")
trachea.one <- CreateSeuratObject(
  counts = counts.trachea.one$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
trachea.one <- AddMetaData(trachea.one, "Trachea", col.name = "group")
trachea.one <- AddMetaData(trachea.one, "tracheaone", col.name = "batch")
trachea.one <- RenameCells(trachea.one, add.cell.id = "T1")

counts.primary.one <- ReadCB_h5("~/scRNAseq/Data/Multiome/Ansley and Peggy Proximal to Distal scRNA-seq/6519-AC/Sample_6519-AC-2/Cellbender v03/output_filtered.h5")
primary.one <- CreateSeuratObject(
  counts = counts.primary.one$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
primary.one <- AddMetaData(primary.one, "primary", col.name = "group")
primary.one <- AddMetaData(primary.one, "primaryone", col.name = "batch")
primary.one <- RenameCells(primary.one, add.cell.id = "P1")

counts.small.one <- ReadCB_h5("~/scRNAseq/Data/Multiome/Ansley and Peggy Proximal to Distal scRNA-seq/6519-AC/Sample_6519-AC-3/Cellbender v03/output_filtered.h5")
small.one <- CreateSeuratObject(
  counts = counts.small.one$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
small.one <- AddMetaData(small.one, "small", col.name = "group")
small.one <- AddMetaData(small.one, "smallone", col.name = "batch")
small.one <- RenameCells(small.one, add.cell.id = "S1")

counts.distal.one <- ReadCB_h5("~/scRNAseq/Data/Multiome/Ansley and Peggy Proximal to Distal scRNA-seq/6519-AC/Sample_6519-AC-4/Cellbender v03/output_filtered.h5")
distal.one <- CreateSeuratObject(
  counts = counts.distal.one$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
distal.one <- AddMetaData(distal.one, "distal", col.name = "group")
distal.one <- AddMetaData(distal.one, "distalone", col.name = "batch")
distal.one <- RenameCells(distal.one, add.cell.id = "D1")

counts.trachea.two <- ReadCB_h5("~/scRNAseq/Data/Multiome/Ansley and Peggy Proximal to Distal scRNA-seq/6534-AC/Sample_6534-AC-1/Cellbender v03/output_filtered.h5")
trachea.two <- CreateSeuratObject(
  counts = counts.trachea.two$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
trachea.two <- AddMetaData(trachea.two, "Trachea", col.name = "group")
trachea.two <- AddMetaData(trachea.two, "tracheatwo", col.name = "batch")
trachea.two <- RenameCells(trachea.two, add.cell.id = "T2")

counts.primary.two <-ReadCB_h5("your/path")
primary.two <- CreateSeuratObject(
  counts = counts.primary.two$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
primary.two <- AddMetaData(primary.two, "primary", col.name = "group")
primary.two <- AddMetaData(primary.two, "primarytwo", col.name = "batch")
primary.two <- RenameCells(primary.two, add.cell.id = "P2")

counts.small.two <- ReadCB_h5("your/path")
small.two <- CreateSeuratObject(
  counts = counts.small.two$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
small.two <- AddMetaData(small.two, "small", col.name = "group")
small.two <- AddMetaData(small.two, "smalltwo", col.name = "batch")
small.two <- RenameCells(small.two, add.cell.id = "S2")

counts.distal.two <-ReadCB_h5("your/path")
distal.two <- CreateSeuratObject(
  counts = counts.distal.two$`Gene Expression`,
  assay = "RNA", min.cells = 3, min.features = 200
)
distal.two <- AddMetaData(distal.two, "distal", col.name = "group")
distal.two <- AddMetaData(distal.two, "distaltwo", col.name = "batch")
distal.two <- RenameCells(distal.two, add.cell.id = "D2")

counts.trachea.three <- ReadCB_h5("your/path")
trachea.three <- CreateSeuratObject(
  counts = counts.trachea.three,
  assay = "RNA", min.cells = 3, min.features = 200
)
trachea.three <- AddMetaData(trachea.three, "Trachea", col.name = "group")
trachea.three <- AddMetaData(trachea.three, "tracheathree", col.name = "batch")
trachea.three <- RenameCells(trachea.three, add.cell.id = "T3")

counts.primary.three <- ReadCB_h5("your/path")
primary.three <- CreateSeuratObject(
  counts = counts.primary.three,
  assay = "RNA", min.cells = 3, min.features = 200
)
primary.three <- AddMetaData(primary.three, "primary", col.name = "group")
primary.three <- AddMetaData(primary.three, "primarythree", col.name = "batch")
primary.three <- RenameCells(primary.three, add.cell.id = "P3")

counts.small.three <- ReadCB_h5("your/path")
small.three <- CreateSeuratObject(
  counts = counts.small.three,
  assay = "RNA", min.cells = 3, min.features = 200
)
small.three <- AddMetaData(small.three, "small", col.name = "group")
small.three <- AddMetaData(small.three, "smallthree", col.name = "batch")
small.three <- RenameCells(small.three, add.cell.id = "S3")

counts.distal.three <- ReadCB_h5("your/path")
distal.three <- CreateSeuratObject(
  counts = counts.distal.three,
  assay = "RNA", min.cells = 3, min.features = 200
)
distal.three <- AddMetaData(distal.three, "distal", col.name = "group")
distal.three <- AddMetaData(distal.three, "distalthree", col.name = "batch")
distal.three <- RenameCells(distal.three, add.cell.id = "D3")

counts.trachea.four <- ReadCB_h5("your/path")
trachea.four <- CreateSeuratObject(
  counts = counts.trachea.four,
  assay = "RNA", min.cells = 3, min.features = 200
)
trachea.four <- AddMetaData(trachea.four, "Trachea", col.name = "group")
trachea.four <- AddMetaData(trachea.four, "tracheafour", col.name = "batch")
trachea.four <- RenameCells(trachea.four, add.cell.id = "T4")

counts.primary.four <- ReadCB_h5("your/path")
primary.four <- CreateSeuratObject(
  counts = counts.primary.four,
  assay = "RNA", min.cells = 3, min.features = 200
)
primary.four <- AddMetaData(primary.four, "primary", col.name = "group")
primary.four <- AddMetaData(primary.four, "primaryfour", col.name = "batch")
primary.four <- RenameCells(primary.four, add.cell.id = "P4")

counts.small.four <- ReadCB_h5("your/path")
small.four <- CreateSeuratObject(
  counts = counts.small.four,
  assay = "RNA", min.cells = 3, min.features = 200
)
small.four <- AddMetaData(small.four, "small", col.name = "group")
small.four <- AddMetaData(small.four, "smallfour", col.name = "batch")
small.four <- RenameCells(small.four, add.cell.id = "S4")

counts.distal.four <- ReadCB_h5("your/path")
distal.four <- CreateSeuratObject(
  counts = counts.distal.four,
  assay = "RNA", min.cells = 3, min.features = 200
)
distal.four <- AddMetaData(distal.four, "distal", col.name = "group")
distal.four <- AddMetaData(distal.four, "distalfour", col.name = "batch")
distal.four <- RenameCells(distal.four, add.cell.id = "D4")

all.combined <- merge(trachea.three, c(primary.three, small.three, trachea.four, primary.four, small.four, distal.four, primary.one, small.one, distal.one, trachea.two, primary.two, small.two, distal.two, trachea.one)) #distal three did not pass QC

all.combined[["percent.ribo"]] <- PercentageFeatureSet(all.combined, pattern = "^RP[SL][[:digit:]]")
all.combined[["percent.mt"]] <- PercentageFeatureSet(all.combined, pattern = "^MT-")
pdf(file.path("./", paste0("Adult Normal Pre-QC ", ".pdf")), w=11, h=8.5)
VlnPlot(all.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0)
dev.off()
pdf(file.path("./", paste0("Adult Normal Pre-QC By Batch", ".pdf")), w=11, h=8.5)
VlnPlot(all.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), split.by = "batch", ncol = 4, pt.size = 0)
dev.off()

all.combined <- subset(all.combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5 & percent.ribo < 7.5 & nCount_RNA > 1000)
pdf(file.path("./", paste0("Adult Normal Post-QC ", ".pdf")), w=11, h=8.5)
VlnPlot(all.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0)
dev.off()
DefaultAssay(all.combined) <- "RNA"
all.combined <- NormalizeData(all.combined)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.combined <- CellCycleScoring(all.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


all.combined <- SplitObject(all.combined, split.by = "batch")
i = 1
for (i in seq_along(all.combined)) {
  all.combined[[i]] <- NormalizeData(all.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(all.combined)
for (i in seq_along(along.with = all.combined)) {
  all.combined[[i]] <- ScaleData(all.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(all.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
all.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(all.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.combined.integrated <- ScaleData(all.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
all.combined.integrated <- RunPCA(all.combined.integrated)
ElbowPlot(all.combined.integrated, ndims = 50)
all.combined.integrated <- RunUMAP(all.combined.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
all.combined.integrated <- FindNeighbors(all.combined.integrated, dims = 1:30)
all.combined.integrated <- FindClusters(all.combined.integrated, resolution = 0.75)
all.combined.integrated$batch <- factor(all.combined.integrated$batch, levels = c("tracheaone", "tracheatwo", "tracheathree", "tracheafour", "primaryone", "primarytwo", "primarythree", "primaryfour", "smallone", "smalltwo", "smallthree", "smallfour", "distalone", "distaltwo","distalfour"))
all.combined.integrated$group <- factor(all.combined.integrated$group, levels = c("Trachea", "primary", "small", "distal"))
pdf(file.path("./", paste0("RNA Louvain", ".pdf")), w=11, h=8.5)
p2 <- DimPlot(all.combined.integrated,  label = TRUE, pt.size = 0.5)
p2 
dev.off()

all.combined.integrated.minustwentynine <- subset(all.combined.integrated, idents = c(29), invert = TRUE)
all.combined.integrated.minustwentynine <- SplitObject(all.combined.integrated.minustwentynine, split.by = "batch")
i = 1
for (i in seq_along(all.combined.integrated.minustwentynine)) {
  all.combined.integrated.minustwentynine[[i]] <- NormalizeData(all.combined.integrated.minustwentynine[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(all.combined.integrated.minustwentynine)
for (i in seq_along(along.with = all.combined.integrated.minustwentynine)) {
  all.combined.integrated.minustwentynine[[i]] <- ScaleData(all.combined.integrated.minustwentynine[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(all.combined.integrated.minustwentynine, anchor.features = features, reduction = "rpca", dims = 1:30)
all.combined.integrated.minustwentynine <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(all.combined.integrated.minustwentynine) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.combined.integrated.minustwentynine <- ScaleData(all.combined.integrated.minustwentynine, vars.to.regress = c("S.Score", "G2M.Score"))
all.combined.integrated.minustwentynine <- RunPCA(all.combined.integrated.minustwentynine)
ElbowPlot(all.combined.integrated.minustwentynine, ndims = 50)
all.combined.integrated.minustwentynine <- RunUMAP(all.combined.integrated.minustwentynine, dims = 1:24, reduction.name = "umap", return.model = TRUE)
all.combined.integrated.minustwentynine <- FindNeighbors(all.combined.integrated.minustwentynine, dims = 1:24)
all.combined.integrated.minustwentynine <- FindClusters(all.combined.integrated.minustwentynine, resolution = 0.75)

pdf(file.path("./", paste0("RNA Louvain", ".pdf")), w=11, h=8.5)
p2 <- DimPlot(all.combined.integrated.minustwentynine,  label = TRUE, pt.size = 0.5)
p2 
dev.off()

all.combined.integrated.minustwentynine.minuseleven <- subset(all.combined.integrated.minustwentynine, idents = c(11), invert = TRUE)
all.combined.integrated.minustwentynine.minuseleven <- SplitObject(all.combined.integrated.minustwentynine.minuseleven, split.by = "batch")
i = 1
for (i in seq_along(all.combined.integrated.minustwentynine.minuseleven)) {
  all.combined.integrated.minustwentynine.minuseleven[[i]] <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(all.combined.integrated.minustwentynine.minuseleven)
for (i in seq_along(along.with = all.combined.integrated.minustwentynine.minuseleven)) {
  all.combined.integrated.minustwentynine.minuseleven[[i]] <- ScaleData(all.combined.integrated.minustwentynine.minuseleven[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(all.combined.integrated.minustwentynine.minuseleven, anchor.features = features, reduction = "rpca", dims = 1:30)
all.combined.integrated.minustwentynine.minuseleven <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.combined.integrated.minustwentynine.minuseleven <- ScaleData(all.combined.integrated.minustwentynine.minuseleven, vars.to.regress = c("S.Score", "G2M.Score"))
all.combined.integrated.minustwentynine.minuseleven <- RunPCA(all.combined.integrated.minustwentynine.minuseleven)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven, dims = 1:27, reduction.name = "umap", return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven, dims = 1:27)
all.combined.integrated.minustwentynine.minuseleven <- FindClusters(all.combined.integrated.minustwentynine.minuseleven, resolution = 1.3)
all.combined.integrated.minustwentynine.minuseleven$batch <- factor(all.combined.integrated.minustwentynine.minuseleven$batch, levels = c("tracheaone", "tracheatwo", "tracheathree", "tracheafour", "primaryone", "primarytwo", "primarythree", "primaryfour", "smallone", "smalltwo", "smallthree", "smallfour", "distalone", "distaltwo", "distalthree", "distalfour"))
all.combined.integrated.minustwentynine.minuseleven$group <- factor(all.combined.integrated.minustwentynine.minuseleven$group, levels = c("Trachea", "primary", "small", "distal"))

pdf(file.path("./", paste0("RNA Louvain", ".pdf")), w=11, h=8.5)
p2 <- DimPlot(all.combined.integrated.minustwentynine.minuseleven,  label = TRUE, pt.size = 0.5)
p2 
dev.off()


trimmed.features = rownames(all.combined.integrated.minustwentynine.minuseleven)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)


#Epithelium subclustering
all.combined.integrated.minustwentynine.minuseleven.epithelium <- subset(all.combined.integrated.minustwentynine.minuseleven, idents = c(6, 14, 17, 8, 10, 19, 21, 26, 35, 24, 16, 36, 34, 33))
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.epithelium) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.epithelium <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.epithelium, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.epithelium, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.epithelium, reduction = "pca", dims = 1:15, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.epithelium, dims = 1:15)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.epithelium, resolution = 0.6)
pdf(file.path("./", paste0("RNA Epithelial Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.epithelium, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.epithelium) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.epithelium <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.epithelium)


cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.epithelium$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.epithelium, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}




#subset small airway and distal enriched epithelial cell types
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.epithelium, resolution = 0.7)
all.combined.integrated.minustwentynine.minuseleven.epithelium.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.epithelium, "5" = "Goblet", "6" = "Serous", "16" = "Deuterosome", "9" = "Multiciliated 1", "2" = "Multiciliated 2", "14" = "GRP+ PNEC", "15" = "GHRL+/RFX6+ PNEC", "17" = "GAP43+/ALK+ Neuronal-like")
all.combined.integrated.minustwentynine.minuseleven.epithelium.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.epithelium, "5" = "Airway Epithelium", "6" = "Sub-Mucosal Gland", "16" = "Multiciliated", "9" = "Multiciliated", "2" = "Multiciliated", "14" = "Neuroendocrine", "15" = "Neuroendocrine", "17" = "Neuroendocrine")

DimPlot(all.combined.integrated.minustwentynine.minuseleven.epithelium, group.by = c("ident"), label = TRUE, pt.size = 0.5)
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- subset(all.combined.integrated.minustwentynine.minuseleven.epithelium, idents = c(1,7,3,4))
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, reduction = "pca", dims = 1:9, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, dims = 1:9)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.epithelium, resolution = 0.7)
pdf(file.path("./", paste0("RNA seclapstlktip Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}

all.combined.integrated.minustwentynine.minuseleven.seclapstlktip.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip,"6" = "Proliferative Distal Epithelium", "5" = "Secretory" ,"4" = "Bud Tip Progenitor 2", "3" = "Bud Tip Progenitor 1", "2" = "Bud Tip Adjacent", "1" = "Lower Airway Progenitor", "0" = "Lower Airway Progenitor")
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip,"6" = "Distal Epithelium", "5" = "Airway Epithelium" , "4" = "Distal Epithelium", "3" = "Distal Epithelium", "2" = "Distal Epithelium", "1" = "Small Airway Epithelium", "0" = "Small Airway Epithelium")
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, col.name = "annotation_lvl2", Idents(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip.named2))
all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, col.name = "annotation_lvl1", "Epithelium")

all.combined.integrated.minustwentynine.minuseleven.seclapstlktip <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip, col.name = "annotation_lvl3", Idents(all.combined.integrated.minustwentynine.minuseleven.seclapstlktip.named))


#subset basal and proximally enriched subtypes
all.combined.integrated.minustwentynine.minuseleven.basal <- subset(all.combined.integrated.minustwentynine.minuseleven.epithelium, idents = c(6, 12, 0, 8, 11))
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.basal <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.basal, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.basal, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.basal <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.basal, reduction = "pca", dims = 1:7, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.basal <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.basal, dims = 1:7)
all.combined.integrated.minustwentynine.minuseleven.basal <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.basal, resolution = 0.2)  
pdf(file.path("./", paste0("RNA basal Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.basal, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.basal <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.basal)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.basal$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.basal, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}
#remove outlier contaminated cluster 4 - unidentifiable
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- subset(all.combined.integrated.minustwentynine.minuseleven.basal, idents = c(4), invert = TRUE)
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, reduction = "pca", dims = 1:6, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, dims = 1:6)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, resolution = 0.3)

pdf(file.path("./", paste0("RNA basal Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}

#remove mesenchymal contaminated cluster 1
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- subset(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour, idents = c(1), invert = TRUE)
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, reduction = "pca", dims = 1:6, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, dims = 1:6)
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, resolution = 0.55)  
pdf(file.path("./", paste0("RNA basal Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone)


cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}

all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, "5" = "KRT4+ Epibasal", "4" = "KRT13+ Epibasal", "3" = "Proliferative Basal", "2" = "LGR5+ Basal 2", "1" = "Transitional Basal", "0" = "LGR5+ Basal 1")
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone,  "5" = "Epibasal" , "4" = "Epibasal", "3" = "Proliferative Basal", "2" = "Small Airway Basal", "1" = "Basal", "0" = "Basal")
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, col.name = "annotation_lvl2", Idents(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone.named2))
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, col.name = "annotation_lvl1", "Epithelium")
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, col.name = "annotation_lvl3", Idents(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone.named))
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3 <- factor(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3, levels = c("LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal"))

all.combined.integrated.basal <- all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone #final basal cell object
#append all epithelial annotations to epithelial object
all.combined.integrated.minustwentynine.minuseleven.epithelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.epithelium, col.name = "annotation_lvl3", c(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.seclapstlktip$annotation_lvl3, Idents(all.combined.integrated.minustwentynine.minuseleven.epithelium.named)))
all.combined.integrated.minustwentynine.minuseleven.epithelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.epithelium, col.name = "annotation_lvl2", c(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl2, all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl2, all.combined.integrated.minustwentynine.minuseleven.seclapstlktip$annotation_lvl2, Idents(all.combined.integrated.minustwentynine.minuseleven.epithelium.named2)))
all.combined.integrated.minustwentynine.minuseleven.epithelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.epithelium, col.name = "annotation_lvl1", c("Epithelium"))




##mesenchymal subclustering
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- subset(all.combined.integrated.minustwentynine.minuseleven, idents = c(0, 1, 5, 4, 3, 2, 31, 23, 15, 11, 22, 7, 9))
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.mesenchyme) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, reduction = "pca", dims = 1:18, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, dims = 1:18)
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, resolution = 0.85)  
pdf(file.path("./", paste0("RNA Mesenchyme Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.mesenchyme) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.mesenchyme)

trimmed.features = rownames(all.combined.integrated.minustwentynine.minuseleven.mesenchyme)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.mesenchyme$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}

all.combined.integrated.minustwentynine.minuseleven.mesenchyme.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.mesenchyme,  '0' = "RSPO2+ Mesenchyme 1", '1' = "RSPO2+ Mesenchyme 2", '2' = "TWIST2+ Mesenchyme", '3' = "GAS2/IGF1 Circumferential Mesenchyme", '4' = "RSPO2+ Mesenchyme 3", '5' = "SCARA5+ Mesenchyme",  '6' = "Fibroblast", "7" = "RSPO2+ Mesenchyme 4", '8' = "PTCH2+ Mesenchyme 2", '9' = "Proliferative Mesenchyme 1", '10' = "Chondrocyte", '11' = "Myofibroblast", '12' = "Chondrocyte Precursor", '13' = "Chondrocyte Precursor", '14' = "PTCH2+ Mesenchyme 1", '15' = "Mature Chondrocyte", 
                                                                                     '16' = "Smooth Muscle", '17' = "Proliferative Mesenchyme 2")
all.combined.integrated.minustwentynine.minuseleven.mesenchyme.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.mesenchyme,  '0' = "Distal Mesenchyme", '1' = "Distal Mesenchyme", '2' = "Proximal Mesenchyme", '3' = "Chondrocyte Lineage", '4' = "Distal Mesenchyme", '5' = "Small Airway Mesenchyme",  '6' = "Mesenchyme", "7" = "Distal Mesenchyme", '8' = "Small Airway Adventitial Mesenchyme", '9' = "Proliferative Mesenchyme", '10' = "Chondrocyte Lineage", '11' = "Smooth Muscle", '12' = "Chondrocyte Lineage", '13' = "Chondrocyte Lineage", '14' = "Proximal Mesenchyme", '15' = "Chondrocyte Lineage", 
                                                                                      '16' = "Smooth Muscle", '17' = "Proliferative Mesenchyme")
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, col.name = "annotation_lvl2", Idents(all.combined.integrated.minustwentynine.minuseleven.mesenchyme.named2))
all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, col.name = "annotation_lvl1", "Mesenchyme")

all.combined.integrated.minustwentynine.minuseleven.mesenchyme <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.mesenchyme, col.name = "annotation_lvl3", Idents(all.combined.integrated.minustwentynine.minuseleven.mesenchyme.named))


all.combined.integrated.mesenchyme <- all.combined.integrated.minustwentynine.minuseleven.mesenchyme #final mesenchyme object
#Endothelium subclustering (also contains pericytes and vsmcs)
all.combined.integrated.minustwentynine.minuseleven.endothelium <- subset(all.combined.integrated.minustwentynine.minuseleven, idents = c(13, 18, 25, 12))
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.endothelium) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.endothelium <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.endothelium, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.endothelium, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.endothelium <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.endothelium, reduction = "pca", dims = 1:5, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.endothelium <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.endothelium, dims = 1:5)
all.combined.integrated.minustwentynine.minuseleven.endothelium <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.endothelium, resolution = 0.125)  
pdf(file.path("./", paste0("RNA endothelium Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.endothelium, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.endothelium) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.endothelium <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.endothelium)

trimmed.features = rownames(all.combined.integrated.minustwentynine.minuseleven.endothelium)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.endothelium$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.endothelium, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}

all.combined.integrated.minustwentynine.minuseleven.endothelium.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.endothelium,  "4" = "Pericytes", "3" = "Arterial Endothelium", "2" = "Vascular Smooth Muscle", "1" = "Capillary Endothelium", "0" = "Lymphatic Endothelium")
all.combined.integrated.minustwentynine.minuseleven.endothelium.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.endothelium,  "4" = "Mural", "3" = "Endothelial", "2" = "Mural", "1" = "Endothelial", "0" = "Endothelial")
all.combined.integrated.minustwentynine.minuseleven.endothelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.endothelium, col.name = "annotation_lvl2", Idents(all.combined.integrated.minustwentynine.minuseleven.endothelium.named2))
all.combined.integrated.minustwentynine.minuseleven.endothelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.endothelium, col.name = "annotation_lvl1", "Vessels")

all.combined.integrated.minustwentynine.minuseleven.endothelium <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.endothelium, col.name = "annotation_lvl3", Idents(all.combined.integrated.minustwentynine.minuseleven.endothelium.named))
all.combined.integrated.minustwentynine.minuseleven.endothelium$annotation_lvl3 <- factor(all.combined.integrated.minustwentynine.minuseleven.endothelium$annotation_lvl3, levels = c("Vascular Smooth Muscle", "Pericytes", "Arterial Endothelium",  "Capillary Endothelium", "Lymphatic Endothelium"))

#Immune subclustering
all.combined.integrated.minustwentynine.minuseleven.immune <- subset(all.combined.integrated.minustwentynine.minuseleven, idents = c(20, 38, 29, 37, 32, 27))
all.combined.integrated.minustwentynine.minuseleven.immune <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.immune, col.name = "annotation_lvl1", "Immune")
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.immune) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.immune <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.immune, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.immune, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.immune <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.immune, reduction = "pca", dims = 1:10, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.immune <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.immune, dims = 1:10)
all.combined.integrated.minustwentynine.minuseleven.immune <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.immune, resolution = 0.15)  
pdf(file.path("./", paste0("RNA immune Louvain", ".pdf")), w=11, h=8.5)
DimPlot(all.combined.integrated.minustwentynine.minuseleven.immune, group.by = c("ident"), label = TRUE, pt.size = 0.5)
dev.off()

DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.immune) <- "RNA"
all.combined.integrated.minustwentynine.minuseleven.immune <- NormalizeData(all.combined.integrated.minustwentynine.minuseleven.immune)

trimmed.features = rownames(all.combined.integrated.minustwentynine.minuseleven.immune)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)

cluster.ids <- levels(all.combined.integrated.minustwentynine.minuseleven.immune$seurat_clusters)
for (i in seq_along(cluster.ids)) {
  markers <- FindMarkers(all.combined.integrated.minustwentynine.minuseleven.immune, ident.1 = cluster.ids[i], min.pct = 0.25, features = trimmed.features)
  write.csv(markers, paste0(cluster.ids[i], "MarkersNoMTorRP.csv"))
}
all.combined.integrated.minustwentynine.minuseleven.immune.named <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.immune, "6" = "Basophils", "5" = "Neutrophils", "4" = "DC1", "3" = "DC2", "2" = "B-Cells", "1" = "Monocytes", "0" = "T-Cells")
all.combined.integrated.minustwentynine.minuseleven.immune.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven.immune, "6" = "Myeloid", "5" = "Myeloid", "4" = "Myeloid", "3" = "Myeloid", "2" = "Lymphoid", "1" = "Myeloid", "0" = "Lymphoid")
all.combined.integrated.minustwentynine.minuseleven.immune <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.immune, col.name = "annotation_lvl2", Idents(all.combined.integrated.minustwentynine.minuseleven.immune.named2))

all.combined.integrated.minustwentynine.minuseleven.immune <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven.immune, col.name = "annotation_lvl3", Idents(all.combined.integrated.minustwentynine.minuseleven.immune.named))
all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl3 <- factor(all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl3, levels = c("Basophils", "Neutrophils", "Monocytes",  "DC2", "DC1", "B-Cells", "T-Cells"))
immune.dotplot.markers <- c("SLC18A2", "MAST4", "MS4A2", "CSF3R",  "MXD1", "ARHGAP26","MRC1", "CD163", "HLA-DRA", "CD74", "CD86",  "GPAT3","CLEC10A", "IRF8", "ZNF366",  "CLNK","CLEC9A", "BANK1", "MS4A1","CD19", "CD247", "CD96", "CD3G")

#transfer annotations from subclustering onto the main object
all.combined.integrated.minustwentynine.minuseleven.named2 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven, "28" = "Neural Crest", "30" = "Sub-Mucosal Gland")
all.combined.integrated.minustwentynine.minuseleven.named3 <- RenameIdents(all.combined.integrated.minustwentynine.minuseleven, "28" = "Neural Crest", "30" = "Epithelium")

all.combined.integrated.minustwentynine.minuseleven <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven, col.name = "annotation_lvl3", c(all.combined.integrated.minustwentynine.minuseleven.mesenchyme$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.epithelium$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.endothelium$annotation_lvl3, Idents(all.combined.integrated.minustwentynine.minuseleven.named)))
all.combined.integrated.minustwentynine.minuseleven <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven, col.name = "annotation_lvl2", c(all.combined.integrated.minustwentynine.minuseleven.mesenchyme$annotation_lvl2, all.combined.integrated.minustwentynine.minuseleven.epithelium$annotation_lvl2, all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl2, all.combined.integrated.minustwentynine.minuseleven.endothelium$annotation_lvl2, Idents(all.combined.integrated.minustwentynine.minuseleven.named2)))
all.combined.integrated.minustwentynine.minuseleven <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven, col.name = "annotation_lvl1", c(all.combined.integrated.minustwentynine.minuseleven.mesenchyme$annotation_lvl1, all.combined.integrated.minustwentynine.minuseleven.epithelium$annotation_lvl1, all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl1, all.combined.integrated.minustwentynine.minuseleven.endothelium$annotation_lvl1, Idents(all.combined.integrated.minustwentynine.minuseleven.named3)))
all.combined.integrated.minustwentynine.minuseleven <- AddMetaData(all.combined.integrated.minustwentynine.minuseleven, col.name = "annotation_lvl4", c(all.combined.integrated.minustwentynine.minuseleven.immune$annotation_lvl4))


#remove unannotated cells from the main object - these were removed at various points in the subclustering and are now removed from the full object
Idents(all.combined.integrated.minustwentynine.minuseleven) <- "annotation_lvl3"
all.combined.integrated.minustwentynine.minuselevenremoveunlabelled <- subset(all.combined.integrated.minustwentynine.minuseleven, idents = c(0, 8, 10, 11, 13), invert = TRUE)

all.combined.integrated <- all.combined.integrated.minustwentynine.minuselevenremoveunlabelled #final full sized object
#remove GAP43+ Neuronal Cluster from epithelium and recluster (tricky - clustered with the neuroendocrine cells so came along in epithelial subclustering)
Idents(all.combined.integrated.minustwentynine.minuseleven.epithelium) <- "annotation_lvl3"
all.combined.integrated.minustwentynine.minuseleven.epithelium <- subset(all.combined.integrated.minustwentynine.minuseleven.epithelium, idents = "GAP43+/ALK+ Neuronal-like", invert = TRUE) #65 cells
DefaultAssay(all.combined.integrated.minustwentynine.minuseleven.epithelium) <- "integrated"
all.combined.integrated.minustwentynine.minuseleven.epithelium <- RunPCA(all.combined.integrated.minustwentynine.minuseleven.epithelium, verbose = FALSE)
ElbowPlot(all.combined.integrated.minustwentynine.minuseleven.epithelium, ndims = 50)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- RunUMAP(all.combined.integrated.minustwentynine.minuseleven.epithelium, reduction = "pca", dims = 1:16, return.model = TRUE)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindNeighbors(all.combined.integrated.minustwentynine.minuseleven.epithelium, dims = 1:16)
all.combined.integrated.minustwentynine.minuseleven.epithelium <- FindClusters(all.combined.integrated.minustwentynine.minuseleven.epithelium, resolution = 0.6)  

all.combined.integrated.epithelium <- all.combined.integrated.minustwentynine.minuseleven.epithelium #final epithelial object

#color coding for annotation level 3
identities <- c("Lymphatic Endothelium", "Capillary Endothelium", "Arterial Endothelium", "Pericytes", "Vascular Smooth Muscle", "Neural Crest", "Myoepithelial", "T-Cells", "Basophils", "B-Cells", "DC1", "DC2", "Neutrophils", "Monocytes", "Multiciliated 1", "Multiciliated 2", "Deuterosome", "GAP43+/ALK+ Neuronal-like", "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "Secretory", "Serous", "Lower Airway Progenitor", "Goblet", "Proliferative Distal Epithelium", "Bud Tip Adjacent", "Bud Tip Progenitor 1", "Bud Tip Progenitor 2", "KRT4+ Epibasal", "Proliferative Basal", "LGR5+ Basal 2", "LGR5+ Basal 1", "Transitional Basal", "KRT13+ Epibasal", "RSPO2+ Mesenchyme 4", "RSPO2+ Mesenchyme 3", "RSPO2+ Mesenchyme 2", "RSPO2+ Mesenchyme 1", "Proliferative Mesenchyme 1", "Proliferative Mesenchyme 2", "Smooth Muscle", "Myofibroblast", "TWIST2+ Mesenchyme", "Fibroblast", "GAS2/IGF1 Circumferential Mesenchyme", "SCARA5+ Mesenchyme", "PTCH2+ Mesenchyme 2", "PTCH2+ Mesenchyme 1", "Chondrocyte Precursor", "Chondrocyte", "Mature Chondrocyte")
colors <- c("#0b3911", "#9ee5a4", "#0b2d6c", "#9f04fc", "#bcc5eb", "#399283", "#399283", "#4a10a1", "#84ee15", "#72182e", "#41d8f4", "#eb1138", "#0df38f", "#e33ab9", "#0ca82e", "#2d68c7", "#e3a0fa", "#799d10", "#6961f9", "#dfd945", "#9f04fc", "#f8ba7c", "#9f6c3b", "#f67a59", "#6e7888", "#ad58a2", "#462a09", "#0b3911", "#9ee5a4", "#0b2d6c", "#bcc5eb", "#4a10a1", "#84ee15", "#72182e", "#41d8f4", "#f8ba7c", "#0df38f", "#e33ab9", "#0ca82e", "#2d68c7", "#e3a0fa", "#799d10", "#6961f9", "#462a09", "#eb1138", "#9f6c3b", "#f67a59", "#6e7888", "#ad58a2", "#9f04fc", "#399283")
color_map <- setNames(colors, identities)



#Fig 1A:
DimPlot(all.combined.integrated,  group.by = "annotation_lvl3", label = TRUE, pt.size = 0.5, cols = color_map)
#Fig 1C:
DimPlot(all.combined.integrated,  split.by = "group", label = TRUE, pt.size = 0.5, cols = color_map)

#Fig 2A
DimPlot(all.combined.integrated.epithelium,  group.by = "annotation_lvl3", label = TRUE, pt.size = 0.5, cols = color_map)
#Fig 2B
all.combined.integrated.epithelium$annotation_lvl3 <- factor(all.combined.integrated.epithelium$annotation_lvl3, levels = c("Bud Tip Progenitor 2","Bud Tip Progenitor 1",  "Bud Tip Adjacent",  "Lower Airway Progenitor", "Secretory", "Goblet", "Serous", "Multiciliated 1", "Multiciliated 2", "Deuterosome",  "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "Myoepithelial", "LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal","Proliferative Distal Epithelium"))
epithelium.dotplot.markers <- c( "DMBT1", "ABCA3",  "C3", "SFTPC","TESC", "ETV5","PIK3C2G", "SOX9", "CPM", "EMP2", "SFTA3", "LAMA3", "CLDN18","SFTPB", "SCGB3A2",  "CFTR", "SCGB1A1", "SCGB3A1", "MUC5B", "LTF", "FOXJ1", "CDC20B", "MUC16", "C6", "ASCL1", "CHGA", "GRP", "NEUROD1", "RFX6", "GHRL", "ACTA2", "LGR6", "TP63","EGFR","KRT17","IL33","ADAMTS3","LGR5", "GABRB2","MMP10","KISS1", "KRT14", "KRT15", "KRT5",  "KRT6A","AQP3", "KRT13", "KRT4","MUC4","TOP2A", "MKI67", "PCNA")
DefaultAssay(all.combined.integrated.epithelium) <- "RNA"
all.combined.integrated.epithelium <- NormalizeData(all.combined.integrated.epithelium)
pdf(file.path("./", paste0("RNA Epithelium dotplot", ".pdf")), w=24, h=15)
Dotplot_Zhiwei_Version(all.combined.integrated.epithelium, epithelium.dotplot.markers, "annotation_lvl3") + theme(axis.text.x = element_text(face = "italic", size = 16))
dev.off()
#Fig 2C
DimPlot(all.combined.integrated.epithelium,  split.by = "group", label = TRUE, pt.size = 0.5, cols = color_map)
write.csv(table(all.combined.integrated.epithelium$annotation_lvl3, all.combined.integrated.epithelium$group))

#Fig 2F
DimPlot(all.combined.integrated.basal,  group.by = "annotation_lvl3", label = TRUE, pt.size = 0.5) 
#featureplots
basal.goi <- c("TP63" ,"KRT5", "KRT13", "KRT17",  "LGR5")
for (i in seq_along(basal.goi)) {
  
  print(FeaturePlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, features = basal.goi[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
  
}
#dotplot
all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3 <- factor(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3, levels = c("LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal"))
basal.dotplot.markers <- c("TP63","EGFR","KRT17","IL33","ADAMTS3","LGR5", "GABRB2","MMP10","KISS1", "KRT14", "KRT15", "KRT5",  "KRT6A","AQP3", "KRT13", "KRT4","MUC4", "MUC16", "FGFR1", "TOP2A", "NGFR" )
Dotplot_Zhiwei_Version(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, basal.dotplot.markers, "annotation_lvl3")

#data for stacked bar
write.csv(table(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$annotation_lvl3, all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone$group), "AnnotationLvl3byRegionComposition.csv")

#Fig 3A
DimPlot(all.combined.integrated.mesenchyme,  group.by = "annotation_lvl3", label = TRUE, pt.size = 0.5, cols = color_map)

#Fig 3B
all.combined.integrated.mesenchyme$annotation_lvl3 <- factor(all.combined.integrated.mesenchyme$annotation_lvl3, levels = c("Smooth Muscle", "Myofibroblast", "RSPO2+ Mesenchyme 1",  "RSPO2+ Mesenchyme 2", "RSPO2+ Mesenchyme 3", "RSPO2+ Mesenchyme 4",  "SCARA5+ Mesenchyme", "PTCH2+ Mesenchyme 2", "PTCH2+ Mesenchyme 1","Fibroblast", "TWIST2+ Mesenchyme", "GAS2/IGF1 Circumferential Mesenchyme", "Chondrocyte Precursor", "Chondrocyte", "Mature Chondrocyte", "Proliferative Mesenchyme 1", "Proliferative Mesenchyme 2"))
mesenchyme.dotplot.markers <- c("ACTA2", "MYH11", "TAGLN", "MYOCD", "ACTG2","HHIP","SLC4A4", "PAG1", "LEF1","NEUROD1","PDGFRA", "PIEZO2", "RSPO2","FGFR4",  "PI15", "SCARA5", "PID1", "PTCH2", "JAG1", "KCND2",  "EBF2","AHR", "NAV3", "CCDC102B", "MME", "TWIST2", "FGF13", "SPARC", "ARHGAP26", "COL24A1", "TRPS1","COL11A1", "CNMD", "ACAN", "COL9A1", "EPYC", "TOP2A" )
Dotplot_Zhiwei_Version(all.combined.integrated.mesenchyme, mesenchyme.dotplot.markers, "annotation_lvl3")

#Fig 3C
DimPlot(all.combined.integrated.mesenchyme,  split.by = "group", label = TRUE, pt.size = 0.5, cols = color_map)
#data for stacked bar
write.csv(table(all.combined.integrated.mesenchyme$annotation_lvl3, all.combined.integrated.mesenchyme$group), "AnnotationLvl3byRegionComposition.csv")

#Fig 3D
regional.mes.goi <- c("COL11A1", "PI15", "RSPO2")
for (i in seq_along(basal.goi)) {
  
  print(FeaturePlot(all.combined.integrated.minustwentynine.minuseleven.basal.minusfour.minusone, features = basal.goi[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
  
}


#Fig 4A 
color_code <- data.frame(
  cluster_name = unlist(str_split("Multiciliated 1
Multiciliated 2
Deuterosome
GAP43+/ALK+ Neuronal-like
GRP+ PNEC
GHRL+/RFX6+ PNEC
Secretory
Serous
Lower Airway Progenitor
Goblet
Proliferative Distal Epithelium
Bud Tip Adjacent
Bud Tip Progenitor 1
Bud Tip Progenitor 2
KRT4+ Epibasal
Proliferative Basal
LGR5+ Basal 2
LGR5+ Basal 1
Transitional Basal
KRT13+ Epibasal
RSPO2+ Mesenchyme 4
RSPO2+ Mesenchyme 3
RSPO2+ Mesenchyme 2
RSPO2+ Mesenchyme 1
Proliferative Mesenchyme 1
Proliferative Mesenchyme 2
Smooth Muscle
Myofibroblast
TWIST2+ Mesenchyme
Fibroblast
GAS2/IGF1 Circumferential Mesenchyme
SCARA5+ Mesenchyme
PTCH2+ Mesenchyme 2
PTCH2+ Mesenchyme 1 
Chondrocyte Precursor
Chondrocyte
Mature Chondrocyte
Myoepithelial
PTCH2+ Mesenchyme 1", pattern = "\n")),
  color = unlist(str_split("#0ca82e
#2d68c7
#e3a0fa
#799d10
#6961f9
#dfd945
#9f04fc
#f8ba7c
#9f6c3b
#f67a59
#6e7888
#ad58a2
#462a09
#0b3911
#9ee5a4
#0b2d6c
#bcc5eb
#4a10a1
#84ee15
#72182e
#41d8f4
#f8ba7c
#0df38f
#e33ab9
#0ca82e
#2d68c7
#e3a0fa
#799d10
#6961f9
#462a09
#eb1138
#9f6c3b
#f67a59
#6e7888
#ad58a2
#9f04fc
#399283
#399283
#6e7888", pattern = "\n"))
)

epi.mesn.sub <- subset(all.combined.integrated, subset = annotation_lvl1 %in% c("Epithelium", "Mesenchyme"))
epi.mesn.sub.final <- subset(epi.mesn.sub, subset = batch %in% c("tracheaone", "tracheatwo", "tracheathree", "tracheafour", "primaryone", "primarytwo", "primarythree", "primaryfour"))

# Recluster and rerun umap
DefaultAssay(epi.mesn.sub.final) <- "RNA"

set.seed(888)

epi.mesn.sub.final <- NormalizeData(epi.mesn.sub.final, normalization.method = "LogNormalize", scale.factor =10000)
epi.mesn.sub.final <- RunUMAP(epi.mesn.sub.final, dims = 1:30)
epi.mesn.sub.final <- FindNeighbors(epi.mesn.sub.final, dims = 1:30)
epi.mesn.sub.final <- FindClusters(epi.mesn.sub.final, graph.name = "integrated_snn", resolution = 0.5, algorithm = 4)

expected_clusters <- unlist(str_split("Multiciliated, Neuroendocrine, Sub-Mucosal Gland, Airway Epithelium, Small Airway Epithelium, Distal Epithelium, Proliferative Basal, Small Airway Basal, Basal, Epibasal, Distal Mesenchyme, Mesenchyme, Proximal Mesenchyme, Small Airway Mesenchyme, Small Airway Adventitial Mesenchyme, Proliferative Mesenchyme, Smooth Muscle, Chondrocyte Lineage", pattern = ", "))

epi.mesn.sub.final$annotation_lvl2 <- factor(epi.mesn.sub.final$annotation_lvl2, levels = expected_clusters)

epi.mesn.sub.final <- SetIdent(epi.mesn.sub.final, value = epi.mesn.sub.final$annotation_lvl2)

all.combined.integrated <- colnames(subset(all.combined.integrated, subset = annotation_lvl3 == "GAP43+/ALK+ Neuronal-like"))

epi.mesn.sub.final <- subset(epi.mesn.sub.final, cells = setdiff(colnames(epi.mesn.sub.final), GAP43_pos_cells))

epi.mesn.sub.final <- RunUMAP(epi.mesn.sub.final, dims = 1:30, return.model = TRUE, reduction.name = "umap_new")

annotation_lvl3_new <- intersect(unique(epi.mesn.sub.final$annotation_lvl3), color_code$cluster_name)

color_code <- color_code %>% 
  dplyr::filter(cluster_name %in% annotation_lvl3_new)

epi.mesn.sub.final$annotation_lvl3 <- factor(epi.mesn.sub.final$annotation_lvl3, levels = color_code$cluster_name)

#Fig 4A 
DimPlot(epi.mesn.sub.final, group.by = "annotation_lvl3", cols = color_code$color, label = T, repel = T) + NoLegend()


#Fig 4B 
gene_list_new <- c(unlist(str_split("FZD1, FZD2, FZD3, FZD4, FZD5, FZD6, FZD7, FZD8, FZD9, FZD10, LRP5, LRP6, LGR4, LGR5, LGR6, WNT2B, WNT3A, WNT4, WNT5A, WNT5B, WNT6, WNT7B, WNT9A, RSPO1, RSPO2, RSPO3, WIF1, SFRP1, RNF43, DKK1, ROR1, ROR2", pattern = ", ")))

cluster_order <- unlist(str_split("Multiciliated 1
Multiciliated 2
Deuterosome
GAP43+/ALK+ Neuronal-like
GRP+ PNEC
GHRL+/RFX6+ PNEC
Secretory
Serous
Lower Airway Progenitor
Goblet
Proliferative Distal Epithelium
Bud Tip Adjacent
Bud Tip Progenitor 1
Bud Tip Progenitor 2
KRT4+ Epibasal
Proliferative Basal
LGR5+ Basal 2
LGR5+ Basal 1
Transitional Basal
KRT13+ Epibasal
Myoepithelial
RSPO2+ Mesenchyme 4
RSPO2+ Mesenchyme 3
RSPO2+ Mesenchyme 2
RSPO2+ Mesenchyme 1
Proliferative Mesenchyme 1
Proliferative Mesenchyme 2
Smooth Muscle
Myofibroblast
TWIST2+ Mesenchyme
Fibroblast
GAS2/IGF1 Circumferential Mesenchyme
SCARA5+ Mesenchyme
PTCH2+ Mesenchyme 2
PTCH2+ Mesenchyme 1 
Chondrocyte Precursor
Chondrocyte
Mature Chondrocyte
PTCH2+ Mesenchyme 1", pattern = "\n"))

epi.mesn.sub.final$annotation_lvl3 <- factor(epi.mesn.sub.final$annotation_lvl3, levels = cluster_order)

Idents(epi.mesn.sub.final) <- epi.mesn.sub.final$annotation_lvl3

draw_dotplots <- function(seurat_object, gene_list, idents, gene_list_name, wid, hei){
  
  DefaultAssay(seurat_object) <- "RNA"
  
  pdf(paste0(gene_list_name, "_dotplot.pdf"), width=wid, height=hei)
  print(DotPlot(seurat_object, features = unique(gene_list), dot.scale = 12, dot.min = 0.05, idents = idents) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.8) +theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_colour_gradientn(colours = c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")) + grids(linetype = "dashed"))
  dev.off()
}

#Fig 4B 
draw_dotplots(epi.mesn.sub.final, gene_list_new, 
              setdiff(color_code$cluster_name, c("Myofibroblast", "RSPO2+ Mesenchyme 1", "RSPO2+ Mesenchyme 2", "RSPO2+ Mesenchyme 3", "RSPO2+ Mesenchyme 4")), 
              "epi_mesn_lvl3/new_updated", 22, 15)

#Fig 4E 

small.epi.mesn.sub <- subset(all.combined.integrated, subset = annotation_lvl1_final %in% c("Epithelium", "Mesenchyme") & batch %in% c("smallone", "smalltwo", "smallthree", "smallfour"))

DefaultAssay(small.epi.mesn.sub) <- "RNA"

set.seed(888)

small.epi.mesn.sub <- NormalizeData(small.epi.mesn.sub, normalization.method = "LogNormalize", scale.factor =10000)
small.epi.mesn.sub <- RunUMAP(small.epi.mesn.sub, dims = 1:30)
small.epi.mesn.sub <- FindNeighbors(small.epi.mesn.sub, dims = 1:30)

small.epi.mesn.sub <- subset(small.epi.mesn.sub, subset = annotation_lvl3 %in% setdiff(unique(small.epi.mesn.sub$annotation_lvl3), c("Bud Tip Progenitor 1", "Bud Tip Progenitor 2", "LGR5+ Basal 1", "KRT4+ Epibasal", "Goblet", "Transitional Basal", "Mature Chondrocyte", "KRT13+ Epibasal", "Serous", "Myoepithelial")))

annotation_lvl3 <- data.frame(table(small.epi.mesn.sub$annotation_lvl3))
colnames(annotation_lvl3) <- c("cluster", "cell_count")

annotation_lvl3 <- annotation_lvl3 %>% dplyr::filter(cell_count > 0)
final_clusters <- as.character(annotation_lvl3$cluster)

small.epi.mesn.sub <- subset(small.epi.mesn.sub, subset = annotation_lvl3 %in% final_clusters)

DefaultAssay(small.epi.mesn.sub) <- "RNA"

set.seed(888)

small.epi.mesn.sub <- NormalizeData(small.epi.mesn.sub, normalization.method = "LogNormalize", scale.factor =10000)
small.epi.mesn.sub <- RunUMAP(small.epi.mesn.sub, dims = 1:30)
small.epi.mesn.sub <- FindNeighbors(small.epi.mesn.sub, dims = 1:30)

color_code_sub <- color_code %>% dplyr::filter(
  cluster_name %in% final_clusters
)

small.epi.mesn.sub$annotation_lvl3 = droplevels(small.epi.mesn.sub$annotation_lvl3, exclude = setdiff(levels(small.epi.mesn.sub$annotation_lvl3), final_clusters))

color_code_sub <- sort(color_code_sub, by = cluster_name)

small.epi.mesn.sub$annotation_lvl3 <- factor(small.epi.mesn.sub$annotation_lvl3, levels = color_code_sub$cluster_name)

#Fig 4E 
DimPlot(small.epi.mesn.sub, group.by = "annotation_lvl3", cols = color_code_sub$color) + easy_remove_axes()

#Fig 4F
gene_list_new <- c(unlist(str_split("FZD1, FZD2, FZD3, FZD4, FZD5, FZD6, FZD7, FZD8, FZD9, FZD10, LRP5, LRP6, LGR4, LGR5, LGR6, WNT2B, WNT3A, WNT4, WNT5A, WNT5B, WNT6, WNT7B, WNT9A, RSPO1, RSPO2, RSPO3, WIF1, SFRP1, RNF43, DKK1, ROR1, ROR2", pattern = ", ")))

ordered_clusters <- unlist(str_split("Multiciliated 1
Multiciliated 2
Deuterosome
GRP+ PNEC
GHRL+/RFX6+ PNEC
Secretory
Serous
Lower Airway Progenitor
Goblet
Proliferative Distal Epithelium
Bud Tip Adjacent
Bud Tip Progenitor 1
Bud Tip Progenitor 2
KRT4+ Epibasal
Proliferative Basal
LGR5+ Basal 2
LGR5+ Basal 1
Transitional Basal
KRT13+ Epibasal
RSPO2+ Mesenchyme 4
RSPO2+ Mesenchyme 3
RSPO2+ Mesenchyme 2
RSPO2+ Mesenchyme 1
Proliferative Mesenchyme 1
Proliferative Mesenchyme 2
Smooth Muscle
Myofibroblast
TWIST2+ Mesenchyme
Fibroblast
GAS2/IGF1 Circumferential Mesenchyme
SCARA5+ Mesenchyme
PTCH2+ Mesenchyme 2
PTCH2+ Mesenchyme 1
Chondrocyte Precursor
Chondrocyte
Mature Chondrocyte
Pericytes
Vascular Smooth Muscle", pattern = "\n"))

ordered_clusters <- ordered_clusters[ordered_clusters %in% final_clusters]

small_epi_mesn_sub$annotation_lvl3 <- factor(small_epi_mesn_sub$annotation_lvl3, levels = ordered_clusters)

Idents(small_epi_mesn_sub) <- small_epi_mesn_sub$annotation_lvl3

draw_dotplots <- function(seurat_object, gene_list, gene_list_name, w, h){
  
  DefaultAssay(seurat_object) <- "RNA"
  
  pdf(paste0(gene_list_name, "_dotplot.pdf"), width=w, height=h)
  print(DotPlot(seurat_object, features = unique(gene_list), dot.scale = 12, dot.min = 0.05) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.8) +theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_colour_gradientn(colours = c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")) + grids(linetype = "dashed"))
  dev.off()
}

draw_dotplots(small_epi_mesn_sub, gene_list_new, "small_epi_mesn_plots/WNT_panel", 22, 15)

#Figure 5
#Map In Vitro Cells onto In Vitro Dataset (MapQuery) Transfer Labels from In Vivo to In Vitro Dataset (TransferData)
counts.trachea.invitro.one <- ReadCB_h5("your/path")
trachea.invitro.one <- CreateSeuratObject(
  counts = counts.trachea.invitro.one,
  assay = "RNA", min.cells = 3, min.features = 200
)
trachea.invitro.one <- AddMetaData(trachea.invitro.one, "trachea", col.name = "group")
trachea.invitro.one <- AddMetaData(trachea.invitro.one, "tracheainvitro.one", col.name = "batch")
counts.primary.invitro.one <- ReadCB_h5("your/path")
primary.invitro.one <- CreateSeuratObject(
  counts = counts.primary.invitro.one,
  assay = "RNA", min.cells = 3, min.features = 200
)
primary.invitro.one <- AddMetaData(primary.invitro.one, "primary", col.name = "group")
primary.invitro.one <- AddMetaData(primary.invitro.one, "primaryinvitro.one", col.name = "batch")
counts.small.invitro.one <- ReadCB_h5("your/path")
small.invitro.one <- CreateSeuratObject(
  counts = counts.small.invitro.one,
  assay = "RNA", min.cells = 3, min.features = 200
)
small.invitro.one <- AddMetaData(small.invitro.one, "small", col.name = "group")
small.invitro.one <- AddMetaData(small.invitro.one, "smallinvitro.one", col.name = "batch")
counts.distal.invitro.one <- ReadCB_h5("your/path")
distal.invitro.one <- CreateSeuratObject(
  counts = counts.distal.invitro.one,
  assay = "RNA", min.cells = 3, min.features = 200
)
distal.invitro.one <- AddMetaData(distal.invitro.one, "distal", col.name = "group")
distal.invitro.one <- AddMetaData(distal.invitro.one, "distalinvitro.one", col.name = "batch")
all.invitro.combined <- merge(trachea.invitro.one, c(primary.invitro.one, small.invitro.one, distal.invitro.one))

all.invitro.combined[["percent.mt"]] <- PercentageFeatureSet(all.invitro.combined, pattern = "^MT-")
all.invitro.combined[["percent.ribo"]] <- PercentageFeatureSet(all.invitro.combined, pattern = "^RP[SL][[:digit:]]")

all.invitro.combined <- subset(all.invitro.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)



DefaultAssay(all.invitro.combined) <- "RNA"
all.invitro.combined <- NormalizeData(all.invitro.combined)
all.invitro.combined <- JoinLayers(all.invitro.combined)
all.invitro.combined <- CellCycleScoring(all.invitro.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

all.invitro.combined <- SplitObject(all.invitro.combined, split.by = "batch")

i = 1
for (i in seq_along(all.invitro.combined)) {
  all.invitro.combined[[i]] <- NormalizeData(all.invitro.combined[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(all.invitro.combined, nfeatures = 1000)
for (i in seq_along(along.with = all.invitro.combined)) {
  all.invitro.combined[[i]] <- ScaleData(all.invitro.combined[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(all.invitro.combined, anchor.features = features, reduction = "rpca", dims = 1:30)
all.invitro.combined.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(all.invitro.combined.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.invitro.combined.integrated <- ScaleData(all.invitro.combined.integrated, vars.to.regress = c("S.Score", "G2M.Score"))
all.invitro.combined.integrated <- RunPCA(all.invitro.combined.integrated)
ElbowPlot(all.invitro.combined.integrated, ndims = 50)
all.invitro.combined.integrated <- RunUMAP(all.invitro.combined.integrated, dims = 1:16, reduction.name = "umap", return.model = TRUE)
all.invitro.combined.integrated <- FindNeighbors(all.invitro.combined.integrated, dims = 1:16)
all.invitro.combined.integrated <- FindClusters(all.invitro.combined.integrated, resolution = 0.5)
all.invitro.combined.integrated$group <- factor(all.invitro.combined.integrated$group, levels = c("trachea", "primary", "small", "distal"))


pdf(file.path("./", paste0("RNA Louvain", ".pdf")), w=11, h=8.5)
p2 <- DimPlot(all.invitro.combined.integrated,  label = TRUE, pt.size = 0.5)
p2 
dev.off()

#remove mesenchymal cluster 11
all.invitro.combined.integrated.epithelium <- subset(all.invitro.combined.integrated, idents = 11, invert = TRUE)
all.invitro.combined.integrated.epithelium <- JoinLayers(all.invitro.combined.integrated.epithelium)
DefaultAssay(all.invitro.combined.integrated.epithelium) <- "RNA"
all.invitro.combined.integrated.epithelium <- SplitObject(all.invitro.combined.integrated.epithelium, split.by = "batch")
i = 1
for (i in seq_along(all.invitro.combined.integrated.epithelium)) {
  all.invitro.combined.integrated.epithelium[[i]] <- NormalizeData(all.invitro.combined.integrated.epithelium[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(all.invitro.combined.integrated.epithelium)
for (i in seq_along(along.with = all.invitro.combined.integrated.epithelium)) {
  all.invitro.combined.integrated.epithelium[[i]] <- ScaleData(all.invitro.combined.integrated.epithelium[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(all.invitro.combined.integrated.epithelium, anchor.features = features, reduction = "rpca", dims = 1:30)
all.invitro.combined.integrated.epithelium <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(all.invitro.combined.integrated.epithelium) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

all.invitro.combined.integrated.epithelium <- ScaleData(all.invitro.combined.integrated.epithelium, vars.to.regress = c("S.Score", "G2M.Score"))
all.invitro.combined.integrated.epithelium <- RunPCA(all.invitro.combined.integrated.epithelium)
ElbowPlot(all.invitro.combined.integrated.epithelium, ndims = 50)
all.invitro.combined.integrated.epithelium <- RunUMAP(all.invitro.combined.integrated.epithelium, dims = 1:20, reduction.name = "umap", return.model = TRUE)
all.invitro.combined.integrated.epithelium <- FindNeighbors(all.invitro.combined.integrated.epithelium, dims = 1:20)
all.invitro.combined.integrated.epithelium <- FindClusters(all.invitro.combined.integrated.epithelium, resolution = 0.5)
all.invitro.combined.integrated.epithelium$group <- factor(all.invitro.combined.integrated.epithelium$group, levels = c("trachea", "primary", "small", "distal"))



pdf(file.path("./", paste0("In Vitro Integrated Louvain", ".pdf")), w=11, h=8.5)
p2 <- DimPlot(all.invitro.combined.integrated.epithelium,  label = TRUE, pt.size = 0.5)
p2 
dev.off()

#Label transfer using Seurat's MapQuery function
DefaultAssay(all.combined.integrated.epithelium) <- "integrated"
DefaultAssay(all.invitro.combined.integrated.epithelium) <- "RNA"
all.invitro.combined.integrated.epithelium <- NormalizeData(all.invitro.combined.integrated.epithelium)
common.features <- intersect(rownames(all.combined.integrated.epithelium), rownames(all.invitro.combined.integrated.epithelium))
length(x = common.features)

all.combined.integrated.epithelium.anchors <- FindTransferAnchors(
  reference = all.combined.integrated.epithelium,
  query = all.invitro.combined.integrated.epithelium,
  reference.assay = "integrated",
  query.assay = "RNA",
  features = common.features,
  reference.reduction = "pca",
  k.filter = 200, 
)


invitro.mapped.invivo.epithelium <- MapQuery(
  anchorset = all.combined.integrated.epithelium.anchors,
  query = all.invitro.combined.integrated.epithelium,
  reference = all.combined.integrated.epithelium,
  refdata = "annotation_lvl3",
  reference.reduction = "pca", 
  reduction.model = "umap",
  
)

#Fig 5B
DimPlot(invitro.mapped.invivo.epithelium, group.by = 'predicted.id', pt.size = 0.5, label = TRUE)
basal.goi.2 <- c("TP63" ,"LGR5", "KRT13", "KRT4")
for (i in seq_along(basal.goi)) {
  print(FeaturePlot(invitro.mapped.invivo.epithelium, features = basal.goi.2[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
  }

#Fig 5C
all.combined.integrated.epithelial$id <- 'reference'
invitro.mapped.invivo.epithelium$id <- 'query'
all.combined.integrated.epithelial$group <- "NA"
refquery <- merge(all.combined.integrated.epithelial, invitro.mapped.invivo.epithelium)
refquery[["umap"]] <- merge(all.combined.integrated.epithelial[["umap"]], invitro.mapped.invivo.epithelium[["ref.umap"]])
refquery$region <- factor(refquery$region, levels = c("Trachea", "Bronchi", "NCA", "Distal", NA))
DimPlot(refquery, group.by = 'region', shuffle = TRUE, label = TRUE)



#Figure 6 - Load Adult Proximal Airway Data
setwd("your/path")
YL9940.1.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
setwd("your/path")
YL9940.2.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
setwd("your/path")
YL9940.3.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
setwd("your/path")
YL9940.4.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
setwd("your/path")
YL9853.1.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)
setwd("your/path")
YL9877.5.data <- ReadCB_h5("output_filtered.h5", use.names = TRUE, unique.features = TRUE)

YL9940.1 <- CreateSeuratObject(counts = YL9940.1.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9940.1 <- RenameCells(YL9940.1, add.cell.id = "9940YL1")
YL9940.1[["percent.mt"]] <- PercentageFeatureSet(YL9940.1, pattern = "^MT-")
YL9940.1[["percent.ribo"]] <- PercentageFeatureSet(YL9940.1, pattern = "^RP[SL][[:digit:]]")
YL9940.1 <- AddMetaData(YL9940.1, "9940YL1", "ID")
YL9940.1 <- AddMetaData(YL9940.1, "Trachea", "region")


YL9940.2 <- CreateSeuratObject(counts = YL9940.2.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9940.2 <- RenameCells(YL9940.2, add.cell.id = "9940YL2")
YL9940.2[["percent.mt"]] <- PercentageFeatureSet(YL9940.2, pattern = "^MT-")
YL9940.2[["percent.ribo"]] <- PercentageFeatureSet(YL9940.2, pattern = "^RP[SL][[:digit:]]")
YL9940.2 <- AddMetaData(YL9940.2, "9940YL2", "ID")
YL9940.2 <- AddMetaData(YL9940.2, "Bronchi", "region")


YL9940.3 <- CreateSeuratObject(counts = YL9940.3.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9940.3 <- RenameCells(YL9940.3, add.cell.id = "9940YL3")
YL9940.3[["percent.mt"]] <- PercentageFeatureSet(YL9940.3, pattern = "^MT-")
YL9940.3[["percent.ribo"]] <- PercentageFeatureSet(YL9940.3, pattern = "^RP[SL][[:digit:]]")
YL9940.3 <- AddMetaData(YL9940.3, "9940YL3", "ID")
YL9940.3 <- AddMetaData(YL9940.3, "NCA", "region")


YL9940.4 <- CreateSeuratObject(counts = YL9940.4.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9940.4 <- RenameCells(YL9940.4, add.cell.id = "9940YL4")
YL9940.4[["percent.mt"]] <- PercentageFeatureSet(YL9940.4, pattern = "^MT-")
YL9940.4[["percent.ribo"]] <- PercentageFeatureSet(YL9940.4, pattern = "^RP[SL][[:digit:]]")
YL9940.4 <- AddMetaData(YL9940.4, "9940YL4", "ID")
YL9940.4 <- AddMetaData(YL9940.4, "Distal", "region")

YL9853.1 <- CreateSeuratObject(counts = YL9853.1.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9853.1 <- RenameCells(YL9853.1, add.cell.id = "YL9853YL1")
YL9853.1[["percent.mt"]] <- PercentageFeatureSet(YL9853.1, pattern = "^MT-")
YL9853.1[["percent.ribo"]] <- PercentageFeatureSet(YL9853.1, pattern = "^RP[SL][[:digit:]]")
YL9853.1 <- AddMetaData(YL9853.1, "YL9853YL1", "ID")
YL9853.1 <- AddMetaData(YL9853.1, "Trachea", "region")

YL9877.5 <- CreateSeuratObject(counts = YL9877.5.data$`Gene Expression`,  min.cells = 3, min.features = 200)
YL9877.5 <- RenameCells(YL9877.5, add.cell.id = "9877YL5")
YL9877.5[["percent.mt"]] <- PercentageFeatureSet(YL9877.5, pattern = "^MT-")
YL9877.5[["percent.ribo"]] <- PercentageFeatureSet(YL9877.5, pattern = "^RP[SL][[:digit:]]")
YL9877.5 <- AddMetaData(YL9877.5, "9877YL5", "ID")
YL9877.5 <- AddMetaData(YL9877.5, "NCA", "region")

adult.proxtodistal <- merge(YL9940.1, c(YL9940.2, YL9940.3, YL9853.1, YL9877.5))
adult.proxtodistal <- subset(adult.proxtodistal, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5 & percent.ribo < 7.5 & nCount_RNA > 1000)
adult.proxtodistal <- JoinLayers(adult.proxtodistal)
DefaultAssay(adult.proxtodistal) <- "RNA"
adult.proxtodistal.rpca <- SplitObject(adult.proxtodistal, split.by = "ID")
for (i in seq_along(adult.proxtodistal.rpca)) {
  adult.proxtodistal.rpca[[i]] <- NormalizeData(adult.proxtodistal.rpca[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(adult.proxtodistal.rpca)
for (i in seq_along(along.with = adult.proxtodistal.rpca)) {
  adult.proxtodistal.rpca[[i]] <- ScaleData(adult.proxtodistal.rpca[[i]], features = features) %>% RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(adult.proxtodistal.rpca, anchor.features = features, reduction = "rpca", dims = 1:30)
adult.proxtodistal.rpca.integrated <- IntegrateData(anchors, dims = 1:30)
DefaultAssay(adult.proxtodistal.rpca.integrated) <- "integrated"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
adult.proxtodistal.rpca.integrated <- ScaleData(adult.proxtodistal.rpca.integrated)
adult.proxtodistal.rpca.integrated <- RunPCA(adult.proxtodistal.rpca.integrated)
ElbowPlot(adult.proxtodistal.rpca.integrated, ndims = 50)
adult.proxtodistal.rpca.integrated <- RunUMAP(adult.proxtodistal.rpca.integrated, dims = 1:30, reduction.name = "umap", return.model = TRUE)
adult.proxtodistal.rpca.integrated <- FindNeighbors(adult.proxtodistal.rpca.integrated, dims = 1:30)
adult.proxtodistal.rpca.integrated <- FindClusters(adult.proxtodistal.rpca.integrated, resolution = 0.75)
DefaultAssay(adult.proxtodistal.rpca.integrated) <- "RNA"
adult.proxtodistal.rpca.integrated <- NormalizeData(adult.proxtodistal.rpca.integrated)

#extract adult epithelium
adult.proxtodistal.rpca.integrated.epi <- subset(adult.proxtodistal.rpca.integrated, idents = c(4, 24, 7, 13, 9, 1, 2, 3, 0, 21, 17, 14, 18, 29, 12, 11, 23))
DefaultAssay(adult.proxtodistal.rpca.integrated.epi) <- "integrated"
adult.proxtodistal.rpca.integrated.epi <- ScaleData(adult.proxtodistal.rpca.integrated.epi)
adult.proxtodistal.rpca.integrated.epi <- RunPCA(adult.proxtodistal.rpca.integrated.epi)
ElbowPlot(adult.proxtodistal.rpca.integrated.epi, ndims = 50)
adult.proxtodistal.rpca.integrated.epi <- RunUMAP(adult.proxtodistal.rpca.integrated.epi, dims = 1:16, reduction.name = "umap", return.model = TRUE)
adult.proxtodistal.rpca.integrated.epi <- FindNeighbors(adult.proxtodistal.rpca.integrated.epi, dims = 1:16)
adult.proxtodistal.rpca.integrated.epi <- FindClusters(adult.proxtodistal.rpca.integrated.epi, resolution = 0.6)

DefaultAssay(adult.proxtodistal.rpca.integrated.epi) <- "RNA"
adult.proxtodistal.rpca.integrated.epi <- NormalizeData(adult.proxtodistal.rpca.integrated.epi)

#remove low quality clusters
adult.proxtodistal.rpca.integrated.epi.cleaned <- subset(adult.proxtodistal.rpca.integrated.epi, idents = c(15, 19), invert = TRUE)
DefaultAssay(adult.proxtodistal.rpca.integrated.epi.cleaned) <- "integrated"
adult.proxtodistal.rpca.integrated.epi.cleaned <- ScaleData(adult.proxtodistal.rpca.integrated.epi.cleaned)
adult.proxtodistal.rpca.integrated.epi.cleaned <- RunPCA(adult.proxtodistal.rpca.integrated.epi.cleaned)
ElbowPlot(adult.proxtodistal.rpca.integrated.epi.cleaned, ndims = 50)
adult.proxtodistal.rpca.integrated.epi.cleaned <- RunUMAP(adult.proxtodistal.rpca.integrated.epi.cleaned, dims = 1:13, reduction.name = "umap", return.model = TRUE)
adult.proxtodistal.rpca.integrated.epi.cleaned <- FindNeighbors(adult.proxtodistal.rpca.integrated.epi.cleaned, dims = 1:13)
adult.proxtodistal.rpca.integrated.epi.cleaned <- FindClusters(adult.proxtodistal.rpca.integrated.epi.cleaned, resolution = 0.6)

DefaultAssay(adult.proxtodistal.rpca.integrated.epi.cleaned) <- "RNA"
adult.proxtodistal.rpca.integrated.epi.cleaned <- NormalizeData(adult.proxtodistal.rpca.integrated.epi.cleaned)

#Basal cell subclustering
adult.proxtodistal.rpca.integrated.epi.basal <- subset(adult.proxtodistal.rpca.integrated.epi.cleaned, idents = c(17, 4, 8))
DefaultAssay(adult.proxtodistal.rpca.integrated.epi.basal) <- "integrated"
adult.proxtodistal.rpca.integrated.epi.basal <- ScaleData(adult.proxtodistal.rpca.integrated.epi.basal)
adult.proxtodistal.rpca.integrated.epi.basal <- RunPCA(adult.proxtodistal.rpca.integrated.epi.basal)
ElbowPlot(adult.proxtodistal.rpca.integrated.epi.basal, ndims = 50)
adult.proxtodistal.rpca.integrated.epi.basal <- RunUMAP(adult.proxtodistal.rpca.integrated.epi.basal, dims = 1:10, reduction.name = "umap", return.model = TRUE)
adult.proxtodistal.rpca.integrated.epi.basal <- FindNeighbors(adult.proxtodistal.rpca.integrated.epi.basal, dims = 1:10)
adult.proxtodistal.rpca.integrated.epi.basal <- FindClusters(adult.proxtodistal.rpca.integrated.epi.basal, resolution = 0.4)

DefaultAssay(adult.proxtodistal.rpca.integrated.epi.basal) <- "RNA"
adult.proxtodistal.rpca.integrated.epi.basal <- NormalizeData(adult.proxtodistal.rpca.integrated.epi.basal)
#Fig 6A
DimPlot(adult.proxtodistal.rpca.integrated.epi.basal, group.by = 'seurat_clusters', label = TRUE)
adult.goi <- c("LGR5" ,"LGR6", "KRT13", "TOP2A")
for (i in seq_along(basal.goi)) {
  print(FeaturePlot(adult.proxtodistal.rpca.integrated.epi.basal, features = adult.goi[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
}
#Fig 6B
adult.basal.markers <- c("TP63", "EGFR", "KRT17", "KRT5", "KRT6A", "KRT13", "KRT4", "MUC4", "AQP3", "ADAMTS3", "LGR5", "LGR6", "GABRB2", "MMP10", "KISS1", "KRT14", "MUC16", "FGFR1", "TOP2A", "NGFR")
Dotplot_Zhiwei_Version(adult.proxtodistal.rpca.integrated.epi.basal, adult.basal.markers, "seurat_clusters")

#Fig 6C - Data Analysis
#Map Adult Epithelium onto Fetal Epithelium
DefaultAssay(all.combined.integrated.epithelium) <- "integrated"
DefaultAssay(adult.proxtodistal.rpca.integrated.epi.cleaned) <- "RNA"
adult.proxtodistal.rpca.integrated.epi.cleaned <- JoinLayers(adult.proxtodistal.rpca.integrated.epi.cleaned)
adult.proxtodistal.rpca.integrated.epi.cleaned <- NormalizeData(adult.proxtodistal.rpca.integrated.epi.cleaned)
common.features <- intersect(rownames(all.combined.integrated.epithelium), rownames(adult.proxtodistal.rpca.integrated.epi.cleaned))
length(x = common.features)

all.combined.integrated.epithelium.anchors <- FindTransferAnchors(
  reference = all.combined.integrated.epithelium,
  query = adult.proxtodistal.rpca.integrated.epi.cleaned,
  reference.assay = "integrated",
  query.assay = "RNA",
  features = common.features,
  reference.reduction = "pca",
  k.filter = 200, 
)


adult.mapped.all.combined.integrated.epithelium <- MapQuery(
  anchorset = all.combined.integrated.epithelium.anchors,
  query = adult.proxtodistal.rpca.integrated.epi.cleaned,
  reference = all.combined.integrated.epithelium,
  refdata = "annotation_lvl3",
  reference.reduction = "pca", 
  reduction.model = "umap",
  
)

all.combined.integrated.epithelium$id <- 'reference'
adult.mapped.all.combined.integrated.epithelium$id <- 'query'
refquery <- merge(adult.mapped.all.combined.integrated.epithelium, all.combined.integrated.epithelium)
refquery[["umap"]] <- merge(adult.mapped.all.combined.integrated.epithelium[["ref.umap"]], all.combined.integrated.epithelium[["umap"]])
refquery$region <- factor(refquery$region, levels = c("Trachea", "Bronchi", "NCA", NA))
refquery$predicted.id <- factor(refquery$predicted.id, levels = identities)

#Figure 6C
DimPlot(refquery, group.by = 'predicted.id', shuffle = TRUE, label = TRUE, cols = color_map) + NoAxes() + NoLegend()


#Figure 6D
Idents(fetal.basal.subset) <- "annotation_lvl3"
fetal.basal.subset.allllgrs <- subset(fetal.basal.subset, idents = c("LGR5+ Basal 1", "LGR5+ Basal 2"))



adult.basal.subset <- AddMetaData(adult.basal.subset, col.name = "plot.annotation", "Adult Basal")
fetal.basal.subset.allllgrs <- AddMetaData(fetal.basal.subset.allllgrs, col.name = "plot.annotation", "Fetal LGR5 Basal 1 and 2")

adult.basal.v.fetal.alllgrs <- merge(adult.basal.subset, fetal.basal.subset.allllgrs)
Idents(adult.basal.v.fetal.alllgrs) <- "plot.annotation"
DefaultAssay(adult.basal.v.fetal.alllgrs) <- "RNA"
adult.basal.v.fetal.alllgrs <- NormalizeData(adult.basal.v.fetal.alllgrs)
trimmed.features = rownames(adult.basal.v.fetal.alllgrs)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)

a.markers <- FindMarkers(adult.basal.v.fetal.alllgrs, ident.1 = "Adult Basal", ident.2 = "Fetal LGR5 Basal 1 and 2", 
                         min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
a.markers$Significant <- ifelse(a.markers$p_val_adj < 0.05 & abs(a.markers$avg_log2FC) > 0.5, "Significant", "Not Significant")

a.volcano_data <- a.markers %>%
  rownames_to_column(var = "gene") %>%
  mutate(logP = -log10(p_val_adj))


genes_to_label <- c("LGR5", "LGR6", "KRT5", "TP63", "KRT13", "KISS1", "MMP10")
a.volcano_data$label <- ifelse(a.volcano_data$gene %in% genes_to_label, a.volcano_data$gene, NA)
ggplot(a.volcano_data, aes(x = avg_log2FC, y = logP, color = Significant)) +
  
  # Plot the points
  geom_point() +
  
  # Scale for significant and non-significant points
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  
  # Labels for specific genes with italicized text
  geom_label_repel(aes(label = label),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.size = 0.5,
                   direction = "y",
                   nudge_y = 2,
                   force = 5,
                   max.overlaps = Inf,
                   fontface = "italic") +  # Make gene names italic
  
  
  # Remove plot title and add axis labels
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  
  # Add a black border around the plot
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.0),  # Add black border
    axis.text = element_text(size = 14, face = "bold"),  # Make axis text bold and larger
    axis.title = element_text(size = 16, face = "bold"),  # Make axis labels bold and larger
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust plot margins
  )

# Figure 6F 

adult.proxtodistal.rpca.integrated.epi.mesn.basal <- subset(adult.proxtodistal.rpca.integrated, idents = c(13, 9, 3, 2, 1, 0, 21, 17, 14, 29, 12, 11, 23, 18, 4, 7, 24, 10, 5, 30, 26, 27))
adult.proxtodistal.rpca.integrated.epi.mesn.basal <- subset(adult.proxtodistal.rpca.integrated.epi.mesn.basal, subset = region %in% c("Trachea", "Bronchi"))

adult.proxtodistal.rpca.integrated.epi.mesn.basal$annotation_lvl1 <- case_when(
  adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident %in% c(13, 9, 3, 2, 1, 0, 21, 17, 14, 29, 12, 11, 23, 18) ~ "epithelium",
  adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident %in% c(10, 5, 30, 26, 27) ~ "mesenchyme",
  adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident %in% c(4, 7, 24) ~ "basal",
  .default = "other"
)

adult.proxtodistal.rpca.integrated.epi.mesn.basal <- NormalizeData(adult.proxtodistal.rpca.integrated.epi.mesn.basal, normalization.method = "LogNormalize", scale.factor =10000)
adult.proxtodistal.rpca.integrated.epi.mesn.basal <- RunUMAP(adult.proxtodistal.rpca.integrated.epi.mesn.basal, dims = 1:30)
adult.proxtodistal.rpca.integrated.epi.mesn.basal <- FindNeighbors(adult.proxtodistal.rpca.integrated.epi.mesn.basal, dims = 1:30)
adult.proxtodistal.rpca.integrated.epi.mesn.basal <- FindClusters(adult.proxtodistal.rpca.integrated.epi.mesn.basal, graph.name = "integrated_snn", resolution = 0.5, algorithm = 4)

DimPlot(adult.proxtodistal.rpca.integrated.epi.mesn.basal, group.by = "annotation_lvl1") + easy_remove_axes() | DimPlot(adult.proxtodistal.rpca.integrated.epi.mesn.basal, label = T) + easy_remove_axes()

# Figure 6G 
gene_list_new <- c(unlist(str_split("FZD1, FZD2, FZD3, FZD4, FZD5, FZD6, FZD7, FZD8, FZD9, FZD10, LRP5, LRP6, LGR4, LGR5, LGR6, WNT2B, WNT3A, WNT4, WNT5A, WNT5B, WNT6, WNT7B, WNT9A, RSPO1, RSPO2, RSPO3, WIF1, SFRP1, RNF43, DKK1, ROR1, ROR2", pattern = ", ")))

draw_dotplots <- function(seurat_object, gene_list, gene_list_name, wid, hei){
  
  DefaultAssay(seurat_object) <- "RNA"
  
  pdf(paste0(gene_list_name, "_dotplot.pdf"), width=wid, height=hei)
  print(DotPlot(seurat_object, features = unique(gene_list), dot.scale = 12, dot.min = 0.05) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.8) +theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_colour_gradientn(colours = c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")) + grids(linetype = "dashed"))
  dev.off()
}

order <- rev(levels(adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident))
adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident <- factor(adult.proxtodistal.rpca.integrated.epi.mesn.basal@active.ident, levels = order)

draw_dotplots(adult.proxtodistal.rpca.integrated.epi.mesn.basal, gene_list_new, "adult_GOL_plots/adult_GOL_epi_mesn_basal_sub_wnt_panel", 15, 10)

#Figure S1C
DimPlot(all.combined.integrated, group.by = "annotation_lvl1", pt.size = 0.5)
#Figure S1D
Idents(all.combined.integrated) <- "annotation_lvl1"
DefaultAssay(all.combined.integrated) <- "RNA"
all.combined.integrated <- NormalizeData(all.combined.integrated)
trimmed.features = rownames(all.combined.integrated)
var_regex = '^MT|^RP' # remove MT and RP genes
trimmed.features = grep(var_regex, trimmed.features, invert=T, value=T)
all.combined.integrated.markers <- FindAllMarkers(all.combined.integrated, min.pct = 0.25, logfc.threshold = 0.25, features = trimmed.features)
all.combined.integrated <- ScaleData(all.combined.integrated)
all.combined.integrated.markers %>% group_by(cluster)  %>% slice_min(n = 5, p_val_adj) -> top20.markers
top20.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top20.markers$cluster)
clusteridents <- levels(Idents(all.combined.integrated))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
genes.for.dotplot <- top20.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)
Dotplot_Zhiwei_Version_NoScale(all.combined.integrated, genes.for.dotplot, "annotation_lvl1")

#Fig. S1E
all.combined.integrated$annotation_lvl3 <- factor(all.combined.integrated$annotation_lvl3, levels = c("Bud Tip Progenitor 2","Bud Tip Progenitor 1",  "Bud Tip Adjacent",  "Lower Airway Progenitor", "Secretory", "Goblet", "Serous", "Multiciliated 1", "Multiciliated 2", "Deuterosome",  "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "GAP43+/ALK+ Neuronal-like","LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal","Proliferative Distal Epithelium", "RSPO2+ Mesenchyme 1",  "RSPO2+ Mesenchyme 2", "RSPO2+ Mesenchyme 3", "RSPO2+ Mesenchyme 4",  "SCARA5+ Mesenchyme", "PTCH2+ Mesenchyme 2", "PTCH2+ Mesenchyme 1","Fibroblast", "TWIST2+ Mesenchyme", "GAS2/IGF1 Circumferential Mesenchyme", "Chondrocyte Precursor", "Chondrocyte", "Mature Chondrocyte", "Proliferative Mesenchyme 1", "Proliferative Mesenchyme 2", "Vascular Smooth Muscle", "Pericytes", "Arterial Endothelium",  "Capillary Endothelium", "Lymphatic Endothelium", "Basophils", "Neutrophils", "Monocytes",  "DC2", "DC1", "B-Cells", "T-Cells"))
all.combined.integrated.dotplot.markers <- c("CDH1", "EPCAM", "C3", "DMBT1", "ABCA3", "SFTPC","TESC", "ETV5","PIK3C2G", "SOX9", "CPM", "EMP2", "SFTA3", "LAMA3", "CLDN18","SFTPB", "SCGB3A2",  "CFTR", "SCGB1A1", "SCGB3A1", "MUC5B", "LTF", "FOXJ1", "CDC20B", "MUC16", "C6", "ASCL1", "CHGA", "GRP", "NEUROD1", "RFX6", "GHRL", "GAP43",  "ALK", "DSCAM", "TP63","EGFR","KRT17","IL33","ADAMTS3","LGR5", "GABRB2","MMP10","KISS1", "KRT14", "KRT15", "KRT5","KRT6A", "AQP3", "KRT13", "KRT4","MUC4", "FGFR1", "TOP2A", "VIM", "COL1A1","NGFR",  "ACTA2", "MYH11", "TAGLN", "MYOCD", "ACTG2","HHIP","SLC4A4", "PAG1", "LEF1","NEUROD1","PDGFRA", "PIEZO2", "RSPO2","FGFR4",  "PI15", "SCARA5", "PID1", "PTCH2", "JAG1", "KCND2",  "EBF2","AHR", "NAV3", "CCDC102B", "MME", "TWIST2", "FGF13", "SPARC", "ARHGAP26", "COL24A1", "TRPS1","COL11A1", "CNMD", "ACAN", "COL9A1", "EPYC", "PDGFRB", "MYO1B", "TEX41", "THY1", "PECAM1", "CDH5", "VWF","DKK2", "GJA5",  "CLU", "NOSTRIN", "PLVAP", "CA4", "NRP2", "PROX1",  "PDPN", "PTPRC", "SLC18A2", "MAST4", "MS4A2", "CSF3R",  "MXD1", "ARHGAP26","MRC1", "CD163", "HLA-DRA", "CD74", "CD86",  "GPAT3","CLEC10A", "IRF8", "ZNF366",  "CLNK","CLEC9A", "BANK1", "MS4A1","CD19", "CD247", "CD96", "CD3G")
all.combined.integrated.dotplot.markers <- make.unique(all.combined.integrated.dotplot.markers)
Dotplot_Zhiwei_Version(all.combined.integrated, all.combined.integrated.dotplot.markers, "annotation_lvl3")
#Fig S2
#Load Spatial Data and Transfer Labels from snRNA-seq data onto Spatial Data
Load_Xenium <- function(path){
  data <- ReadXenium(
    data.dir = path,
    type = c("centroids", "segmentations"))
  
  assay <- "Xenium"
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations))
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay)
  
  xenium_obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  xenium_obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium_obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  data$matrix[["Unassigned Codeword"]] <- data$matrix[["Negative Control Codeword"]]
  data$matrix[["Unassigned Codeword"]][data$matrix[["Unassigned Codeword"]] != 0] <- 0
  xenium_obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  
  fov <- "fov"
  xenium_obj[[fov]] <- coords
  
  return(xenium_obj)
}

label_transfer <- function(xenium_object, folder_name, file_name){
  # Find anchors
  set.seed(888)
  xenium_anchors <- FindTransferAnchors(reference = all.combined.integrated, 
                                        query = xenium_object, 
                                        features = xenium_genes, 
                                        normalization.method = "LogNormalize", 
                                        k.filter = 200)
  

  # Transfer labels
  predictions <- TransferData(
    anchorset = xenium_anchors,
    refdata = all.combined.integrated@annotation_lvl3
  )
  
  # Add metadata
  xenium_object <- AddMetaData(object = xenium_object, metadata = predictions)
  
  clusters <- tibble(cluster_name = unique(xenium_object$predicted.id))
  
  color_code <- merge(color_code, clusters, by = "cluster_name")
  
  xenium_object$predicted.id <- factor(xenium_object$predicted.id, levels = color_code$cluster_name)
  
  DefaultAssay(xenium_object) <- "Xenium"
  
  png(paste0(folder_name, "/", file_name, ".png"), width=20, height=12, units="in", res=1024)
  print(ImageDimPlot(xenium_object, group.by = "predicted.id", size = 0.8, cols = color_code$color_list))
  dev.off()
  
  return(xenium_object)
}

cell_ID_group <- function(xenium_object, folder_name, file_name){
  
  df <- data.frame(
    cell_id = names(xenium_object$predicted.id),
    group = xenium_object$predicted.id,
    row.names = NULL
  )
  
  write.csv(df, paste0(folder_name, "/", file_name, ".csv"), row.names = F)
}

xenium.obj.b <- Load_Xenium("11155-TF/output-XETG00077__0040759__11155-TF-3_ROI_B__20240809__163438")
xenium.obj.b <- label_transfer(xenium_ROI_B, "folder", "file")
xenium.obj.b=xenium.obj.b[,unname(which(colSums(GetAssayData(xenium.obj.b))!=0))] #remove cells in which counts log_umi = 0
  
xenium.obj.b <- SCTransform(xenium.obj.b, assay = "Xenium")
xenium.obj.b <- RunPCA(xenium.obj.b, npcs = 30, features = rownames(xenium.obj.b))
  xenium.obj.b <- RunUMAP(xenium.obj.b, dims = 1:30)
  xenium.obj.b <- FindNeighbors(xenium.obj.b, reduction = "pca", dims = 1:30)
  xenium.obj.b <- FindClusters(xenium.obj.b, resolution = 0.3)
  
  
  xenium.obj.d <- Load_Xenium("11155-TF/output-XETG00077__0040756__11155-TF-2_ROI_D__20240809__163438") 
  xenium.obj.d <- label_transfer(xenium_ROI_B, folder, file)
  xenium.obj.d=xenium.obj.d[,unname(which(colSums(GetAssayData(xenium.obj.d))!=0))] #remove cells in which counts log_umi = 0
  
  xenium.obj.d <- SCTransform(xenium.obj.d, assay = "Xenium")
  xenium.obj.d <- RunPCA(xenium.obj.d, npcs = 30, features = rownames(xenium.obj.d))
  xenium.obj.d <- RunUMAP(xenium.obj.d, dims = 1:30)
  xenium.obj.d <- FindNeighbors(xenium.obj.d, reduction = "pca", dims = 1:30)
  xenium.obj.d <- FindClusters(xenium.obj.d, resolution = 0.3)
  
#Fig S2 UMAPS
DimPlot(xenium.obj.b, group.by = "predicted.id", pt.size = 0.5, cols = color_map)
DimPlot(xenium.obj.d, group.by = "predicted.id", pt.size = 0.5, cols = color_map)

#Figure S2 Dotplots
new_order_single <- c("Bud Tip Progenitor 2",
                      "Bud Tip Progenitor 1",
                      "Bud Tip Adjacent",
                      "Lower Airway Progenitor",
                      "Secretory",
                      "Proliferative Distal Epithelium",
                      "Proliferative Basal",
                      "LGR5+ Basal 2",
                      "LGR5+ Basal 1",
                      "Transitional Basal",
                      "KRT13+ Epibasal",
                      "KRT4+ Epibasal",
                      "Myoepithelial",
                      "Serous",
                      "Goblet",
                      "Deuterosome",
                      "Multiciliated 1",
                      "Multiciliated 2",
                      "GRP+ PNEC",
                      "GHRL+/RFX6+ PNEC",
                      "GAP43+/ALK+ Neuronal-like",
                      "Neural Crest",
                      "GAS2/IGF1 Circumferential Mesenchyme",
                      "Chondrocyte Precursor",
                      "Chondrocyte",
                      "Mature Chondrocyte",
                      "Fibroblast",
                      "TWIST2+ Mesenchyme",
                      "PTCH2+ Mesenchyme 1",
                      "PTCH2+ Mesenchyme 2",
                      "SCARA5+ Mesenchyme",
                      "RSPO2+ Mesenchyme 1",
                      "RSPO2+ Mesenchyme 2",
                      "RSPO2+ Mesenchyme 3",
                      "RSPO2+ Mesenchyme 4",
                      "Myofibroblast",
                      "Smooth Muscle",
                      "Vascular Smooth Muscle",
                      "Pericytes",
                      "Proliferative Mesenchyme 1",
                      "Proliferative Mesenchyme 2",
                      "Arterial Endothelium",
                      "Capillary Endothelium",
                      "Lymphatic Endothelium",
                      "Monocytes",
                      "Basophils",
                      "DC1",
                      "DC2",
                      "Neutrophils",
                      "B-Cells",
                      "T-Cells")

xenium.obj.b$predicted.id <- factor(xenium.obj.b$predicted.id, levels = new_order_single)
all.combined.integrated$annotation_lvl3 <- factor(all.combined.integrated$annotation_lvl3, levels = new_order_single)

#Proximal (Top)
DefaultAssay(xenium.obj.b) <- "Xenium"
xenium.obj.b <- NormalizeData(xenium.obj.b)
Idents(xenium.obj.b) <- "predicted.id"
xenium.obj.b.markers <- FindAllMarkers(xenium.obj.b, min.pct = 0.25, logfc.threshold = 0.25)
xenium.obj.b <- ScaleData(xenium.obj.b)
xenium.obj.b.markers %>% group_by(cluster)  %>% top_n(n = 5, wt = avg_log2FC) -> top20.markers
top20.markers$cluster <- gsub(pattern = "/", replacement = "-", x = top20.markers$cluster)
clusteridents <- levels(Idents(xenium.obj.b))
clusteridents <- gsub(pattern = "/", replacement = "-", x = clusteridents)
xenium.obj.b$predicted.id <- factor(xenium.obj.b$predicted.id, levels = sort(unique(xenium.obj.b$predicted.id)))
top20.markers <- top20.markers %>% arrange(cluster)
genes.for.dotplot <- top20.markers$gene
genes.for.dotplot <- make.unique(genes.for.dotplot)

genes.for.dotplot <- genes.for.dotplot[!grepl("\\.", genes.for.dotplot)]


dotplot1 <- Dotplot_Zhiwei_Version_NoScale(all.combined.integrated, genes.for.dotplot,  'annotation_lvl3')
dotplot2 <- Dotplot_Zhiwei_Version_NoScale(xenium.obj.b,  genes.for.dotplot, 'predicted.id')

dotplot1$data$id <- paste("snRNAseq", dotplot1$data$id, sep = "-")
dotplot2$data$id <- paste("Xenium", dotplot2$data$id, sep = "-")
combined_data <- rbind(dotplot1$data, dotplot2$data)
combined_data$id <- factor(combined_data$id, levels = combined.levels)

pdf(file.path("./", paste0("Xenium Proximal (B) Full Top5 Markers per cluster Xenium_All Interleaved", ".pdf")), w=150, h=60)
ggplot(combined_data, aes(x = features.plot, y = id)) + 
  RotatedAxis() +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.4) +  # Add fill = avg.exp.scaled
  scale_size(range = c(1, 20)) +  # Increase the size range of the dots (adjust as needed)
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +  # Use scale_fill_distiller for fill
  guides(size = guide_legend(title = "Percent Expressed", 
                             override.aes = list(shape = 21, colour = "black", fill = "black"))) +
  labs(y = NULL, x = NULL) +
  guides(fill = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")) +
  theme_minimal() +  # Minimal theme for white background
  theme(axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.8)) +
  theme(axis.text.x = element_text(face = "italic", size = 30, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 30)) +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) 
dev.off()


#Distal (Bottom)

dotplot1 <- Dotplot_Zhiwei_Version_NoScale(all.combined.integrated, genes.for.dotplot,  'annotation_lvl3')
dotplot2 <- Dotplot_Zhiwei_Version_NoScale(xenium.obj.d,  genes.for.dotplot, 'predicted.id')

dotplot1$data$id <- paste("snRNAseq", dotplot1$data$id, sep = "-")
dotplot2$data$id <- paste("Xenium", dotplot2$data$id, sep = "-")
combined_data <- rbind(dotplot1$data, dotplot2$data)
combined_data$id <- factor(combined_data$id, levels = combined.levels)

ggplot(combined_data, aes(x = features.plot, y = id)) + 
  RotatedAxis() +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21, colour = "black", stroke = 0.4) +  # Add fill = avg.exp.scaled
  scale_size(range = c(1, 20)) +  # Increase the size range of the dots (adjust as needed)
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +  # Use scale_fill_distiller for fill
  guides(size = guide_legend(title = "Percent Expressed", 
                             override.aes = list(shape = 21, colour = "black", fill = "black"))) +
  labs(y = NULL, x = NULL) +
  guides(fill = guide_colourbar(title = "Average Expression", ticks = TRUE, frame.colour = "black")) +
  theme_minimal() +  # Minimal theme for white background
  theme(axis.line.x.bottom = element_blank(), axis.line.y.left = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.8)) +
  theme(axis.text.x = element_text(face = "italic", size = 30, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 30)) +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())  # Remove minor gridlines

#Fig S3A-B 

PNEC.clusters <- subset(all.combined.integrated.epithelium, idents = c("GRP+ PNEC", "GHRL+/RFX6+ PNEC"))

DefaultAssay(PNEC.clusters) <- "RNA"

set.seed(888)

PNEC.clusters <- NormalizeData(PNEC.clusters, normalization.method = "LogNormalize", scale.factor =10000)
PNEC.clusters <- ScaleData(PNEC.clusters)
PNEC.clusters <- RunUMAP(PNEC.clusters, dims = 1:30)
PNEC.clusters <- FindNeighbors(PNEC.clusters, dims = 1:30)
PNEC.clusters <- FindClusters(PNEC.clusters, graph.name = "integrated_snn",resolution = 0.1, algorithm = 4)

DimPlot(PNEC.clusters)

all.combined.integrated.epithelium.cell.count <- data.frame(table(all.combined.integrated.epithelium$group, all.combined.integrated.epithelium@active.ident)) %>% spread(Var1, Freq)

#Fig S6A
DimPlot(all.invitro.combined.integrated,  group.by = "seurat_clusters", label = TRUE, pt.size = 0.5)
#Fig S6B:
DimPlot(all.invitro.combined.integrated,  split.by = "group", label = TRUE, pt.size = 0.5)
#Fig S6C
write.csv(table(all.invitro.combined.integrated$seurat_clusters, all.invitro.combined.integrated$group))
#Fig S6D
##Find in vitro and in vivo common marker set
all.combined.integrated.epithelial <- in_vivo_new
all.invitro.combined.integrated.epithelium <- transfer_query
DefaultAssay(all.combined.integrated.epithelial) <- "RNA"
DefaultAssay(all.invitro.combined.integrated.epithelium) <- "RNA"
all.combined.integrated.epithelial <- NormalizeData(all.combined.integrated.epithelial)
all.invitro.combined.integrated.epithelium <- NormalizeData(all.invitro.combined.integrated.epithelium)
Idents(all.combined.integrated.epithelial) <- "annotation_lvl3"
Idents(all.invitro.combined.integrated.epithelium) <- "predicted.id"
# Identify markers for each cell type in both datasets
markers_obj1 <- FindAllMarkers(all.combined.integrated.epithelial, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
                               group.by = "annotation_lvl3")

markers_obj2 <- FindAllMarkers(all.invitro.combined.integrated.epithelium, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
                               group.by = "predicted.id")

# Filter top 10 markers for each cluster
top_markers_obj1 <- markers_obj1 %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

top_markers_obj2 <- markers_obj2 %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

# Now, for each cell type (cluster), find common markers between datasets
common_markers <- list()

cell_types <- unique(top_markers_obj1$cluster)

for (cell_type in cell_types) {
  # Get markers for this cell type in both datasets
  markers1 <- top_markers_obj1 %>% filter(cluster == cell_type) %>% pull(gene)
  markers2 <- top_markers_obj2 %>% filter(cluster == cell_type) %>% pull(gene)
  
  # Find common markers
  common <- intersect(markers1, markers2)
  
  # Store the top 3 common markers
  common_markers[[cell_type]] <- head(common, 5)
}

# Display the common markers for each cell type
common_markers
genes_to_plot <- unique(unlist(common_markers))
genes_to_plot
all.combined.integrated.epithelial <- AddMetaData(all.combined.integrated.epithelial, col.name = "plot_annotation", paste0("In Vivo-", Idents(all.combined.integrated.epithelial)))
all.invitro.combined.integrated.epithelium <- AddMetaData(all.invitro.combined.integrated.epithelium, col.name = "plot_annotation", paste0("In Vitro-", Idents(all.invitro.combined.integrated.epithelium)))

all.combined.integrated.epithelial.plot <- merge(all.combined.integrated.epithelial, all.invitro.combined.integrated.epithelium)

all.combined.integrated.epithelial.plot$shared_label <- gsub("In Vivo-|In Vitro-", "", all.combined.integrated.epithelial.plot$plot_annotation)

# Get unique cell type labels and reorder them
unique_labels <- unique(all.combined.integrated.epithelial.plot$shared_label)

# Create a new order for cell types where In Vivo and In Vitro labels are interleaved
new_order <- c()
for (label in unique_labels) {
  in_vivo <- paste("In Vivo-", label, sep = "")
  in_vitro <- paste("In Vitro-", label, sep = "")
  
  # Check if both In Vivo and In Vitro exist for this label
  if (in_vivo %in% all.combined.integrated.epithelial.plot$plot_annotation & in_vitro %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vivo, in_vitro)
  } else if (in_vivo %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vivo)
  } else if (in_vitro %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vitro)
  }
}

# Update the levels of the cell type metadata to match the new order
genes_to_plot <- c("SFTPC",
                   "ROS1",
                   "SLC34A2",
                   "LINC01331",
                   "ROR1",
                   "AC107223.1",
                   "WIF1",
                   "TESC",
                   "PIK3C2G",
                   "MKI67",
                   "KIF23",
                   "TOP2A",
                   "ADAMTSL1",
                   "GDF15",
                   "CFTR",
                   "SCGB3A2",
                   "STEAP4",
                   "CAPN8",
                   "KYNU",
                   "SCGB1A1",
                   "LGR5",
                   "TP63",
                   "COL17A1",
                   "LRP4",
                   "GPC3",
                   "KRT5",
                   "KRT15",
                   "S100A2",
                   "KRT13",
                   "MUC16",
                   "KRT4",
                   "CYP2F1",
                   "CDC20B",
                   "FOXJ1",
                   "ERICH3",
                   "DTHD1",
                   "KIAA2012",
                   "RGS7",
                   "MIAT",
                   "KCNB2", "ASCL1", "GHRL", "GRP")
new_order <- c("In Vivo-Bud Tip Progenitor 2",
               "In Vitro-Bud Tip Progenitor 2",
               "In Vivo-Bud Tip Progenitor 1",
               "In Vitro-Bud Tip Progenitor 1",
               "In Vivo-Bud Tip Adjacent",
               "In Vitro-Bud Tip Adjacent",
               "In Vivo-Proliferative Distal Epithelium",
               "In Vitro-Proliferative Distal Epithelium",
               "In Vivo-Lower Airway Progenitor",
               "In Vitro-Lower Airway Progenitor",
               "In Vivo-Secretory",
               "In Vitro-Secretory",
               "In Vivo-LGR5+ Basal 1",
               "In Vivo-LGR5+ Basal 2",
               "In Vivo-Transitional Basal",
               "In Vivo-KRT13+ Epibasal",
               "In Vitro-KRT13+ Epibasal",
               "In Vivo-KRT4+ Epibasal",
               "In Vitro-KRT4+ Epibasal",
               "In Vivo-Proliferative Basal",
               "In Vitro-Proliferative Basal",
               "In Vivo-Deuterosome",
               "In Vitro-Deuterosome",
               "In Vivo-Multiciliated 1",
               "In Vitro-Multiciliated 1",
               "In Vivo-Multiciliated 2",
               "In Vitro-Multiciliated 2",
               "In Vivo-Serous",
               "In Vitro-Serous",
               "In Vivo-GHRL+/RFX6+ PNEC",
               "In Vivo-GRP+ PNEC",
               "In Vitro-GRP+ PNEC")


all.combined.integrated.epithelial.plot$plot_annotation <- factor(all.combined.integrated.epithelial.plot$plot_annotation, levels = new_order)
Idents(all.combined.integrated.epithelial.plot) <- "plot_annotation"
all.combined.integrated.epithelial.plot <- subset(all.combined.integrated.epithelial.plot, idents = c( "In Vivo-LGR5+ Basal 2", "In Vivo-Goblet"), invert = TRUE)    
all.combined.integrated.epithelial.plot <- NormalizeData(all.combined.integrated.epithelial.plot)
#Plot
Dotplot_Zhiwei_Version_NoScale(all.combined.integrated.epithelial.plot, genes_to_plot, "plot_annotation") + theme(axis.text.x = element_text(face = "italic", size = 16))
#Fig S6E
write.csv(table(all.invitro.combined.integrated.epithelium$predicted.id), "FigS6E.csv")
#Fig S7A
DimPlot(adult.proxtodistal.rpca.integrated.epi.cleaned, group.by = "seurat_clusters", label = TRUE, pt.size = 0.8)
fig7a.goi <- c("TP63" ,"ACTA2")
for (i in seq_along(fig7a.goi)) {
  print(FeaturePlot(adult.proxtodistal.rpca.integrated.epi.cleaned, features = fig7a.goi[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
}
DimPlot(adult.mapped.all.combined.integrated.epithelium, group.by = "predicted.id", label = TRUE, pt.size = 0.8)
write.csv(table(adult.mapped.all.combined.integrated.epithelium$predicted.id), "fig7a.csv")

#Fig S7B
DefaultAssay(all.combined.integrated.epithelial) <- "RNA"
DefaultAssay(adult.proxtodistal.rpca.integrated.epi.cleaned) <- "RNA"
all.combined.integrated.epithelial <- NormalizeData(all.combined.integrated.epithelial)
adult.proxtodistal.rpca.integrated.epi.cleaned <- NormalizeData(adult.proxtodistal.rpca.integrated.epi.cleaned)
adult.proxtodistal.rpca.integrated.epi.cleaned <- AddMetaData(adult.proxtodistal.rpca.integrated.epi.cleaned, col.name = "predicted.id", invitro.mapped.invivo.epithelium$predicted.id)

all.combined.integrated.epithelial$annotation_lvl3 <- factor(all.combined.integrated.epithelial$annotation_lvl3, levels = c("Bud Tip Progenitor 2","Bud Tip Progenitor 1",  "Bud Tip Adjacent",  "Lower Airway Progenitor", "Secretory", "Goblet", "Serous", "Multiciliated 1", "Multiciliated 2", "Deuterosome",  "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "Myoepithelial", "LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal","Proliferative Distal Epithelium"))
adult.proxtodistal.rpca.integrated.epi.cleaned$predicted.id <- factor(adult.proxtodistal.rpca.integrated.epi.cleaned$predicted.id, levels = c("Bud Tip Progenitor 2","Bud Tip Progenitor 1",  "Bud Tip Adjacent",  "Lower Airway Progenitor", "Secretory", "Goblet", "Serous", "Multiciliated 1", "Multiciliated 2", "Deuterosome",  "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "Myoepithelial","LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal","Proliferative Distal Epithelium"))
Idents(all.combined.integrated.epithelial) <- "annotation_lvl3"
Idents(adult.proxtodistal.rpca.integrated.epi.cleaned) <- "predicted.id"
# Identify markers for each cell type in both datasets
markers_obj1 <- FindAllMarkers(all.combined.integrated.epithelial, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
)

markers_obj2 <- FindAllMarkers(adult.proxtodistal.rpca.integrated.epi.cleaned, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
)

# Filter top 10 markers for each cluster
top_markers_obj1 <- markers_obj1 %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

top_markers_obj2 <- markers_obj2 %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)

# Now, for each cell type (cluster), find common markers between datasets
common_markers <- list()

cell_types <- unique(top_markers_obj1$cluster)

for (cell_type in cell_types) {
  # Get markers for this cell type in both datasets
  markers1 <- top_markers_obj1 %>% filter(cluster == cell_type) %>% pull(gene)
  markers2 <- top_markers_obj2 %>% filter(cluster == cell_type) %>% pull(gene)
  
  # Find common markers
  common <- intersect(markers1, markers2)
  
  # Store the top 3 common markers
  common_markers[[cell_type]] <- head(common, 5)
}

# Display the common markers for each cell type
common_markers
genes_to_plot <- unique(unlist(common_markers))
genes_to_plot
all.combined.integrated.epithelial <- AddMetaData(all.combined.integrated.epithelial, col.name = "plot_annotation", paste0("Fetal-", Idents(all.combined.integrated.epithelial)))
adult.proxtodistal.rpca.integrated.epi.cleaned <- AddMetaData(adult.proxtodistal.rpca.integrated.epi.cleaned, col.name = "plot_annotation", paste0("Adult-", Idents(adult.proxtodistal.rpca.integrated.epi.cleaned)))


all.combined.integrated.epithelial.plot <- merge(all.combined.integrated.epithelial, adult.proxtodistal.rpca.integrated.epi.cleaned)

all.combined.integrated.epithelial.plot$shared_label <- gsub("Fetal-|Adult-", "", all.combined.integrated.epithelial.plot$plot_annotation)

# Get unique cell type labels and reorder them
unique_labels <- unique(all.combined.integrated.epithelial.plot$shared_label)
unique_labels_2 <- c("Bud Tip Progenitor 2","Bud Tip Progenitor 1",  "Bud Tip Adjacent",  "Lower Airway Progenitor", "Secretory", "Goblet", "Serous", "Multiciliated 1", "Multiciliated 2", "Deuterosome",  "GRP+ PNEC", "GHRL+/RFX6+ PNEC", "Myoepithelial", "LGR5+ Basal 2","LGR5+ Basal 1",  "Transitional Basal", "KRT13+ Epibasal", "KRT4+ Epibasal", "Proliferative Basal","Proliferative Distal Epithelium")
# Create a new order for cell types where In Vivo and In Vitro labels are interleaved
new_order <- c()
for (label in unique_labels_2) {
  in_vivo <- paste("Fetal-", label, sep = "")
  in_vitro <- paste("Adult-", label, sep = "")
  
  # Check if both In Vivo and In Vitro exist for this label
  if (in_vivo %in% all.combined.integrated.epithelial.plot$plot_annotation & in_vitro %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vivo, in_vitro)
  } else if (in_vivo %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vivo)
  } else if (in_vitro %in% all.combined.integrated.epithelial.plot$plot_annotation) {
    new_order <- c(new_order, in_vitro)
  }
}

# Update the levels of the cell type metadata to match the new order
genes_to_plot <- c("SFTPC", "ROS1", "ABCA3", "AC046195.1", "LAMP3", "ACOXL", "WIF1", "CD36", "GPC5", "ROR1",
                   "KHDRBS2", "LAMA3", "DLC1", "SCGB3A1", "CXCL17", "LCN2", "RIMS1", "SSR4", "KYNU", "BPIFB1",
                   "TGM2", "SLC4A4", "MUC5B", "DACH2", "BPIFB2", "FGF13", "KCNMA1", "LHFPL2", "ITPR2", "CHRM3",
                   "KIAA1324", "CCL28", "CAPS", "TMEM190", "C20orf85", "FOXJ1", "ZEB2", "CFAP157", "ERICH3", "DTHD1",
                   "ZBBX", "LMNTD1", "CDC20B", "CCNO", "KIF24", "FOXN4", "STIL", "ACTA2", "NTRK2", "FGFR1",
                   "LAMA1", "COL4A2", "TP63", "LGR5", "LGR6", "MMP10", "KISS1", "AQP3", "KRT5", "PKP1", "S100A2", "SERPINB13", "ALDH1A3", "LIPH", "EMP1",
                   "DENND2C", "CEACAM6", "TOP2A", "MKI67")

all.combined.integrated.epithelial.plot$plot_annotation <- factor(all.combined.integrated.epithelial.plot$plot_annotation, levels = new_order)
Idents(all.combined.integrated.epithelial.plot) <- "plot_annotation"
all.combined.integrated.epithelial.plot <- subset(all.combined.integrated.epithelial.plot, idents = c( "In Vivo-LGR5+ Basal 2", "In Vivo-Goblet"), invert = TRUE)    
all.combined.integrated.epithelial.plot <- NormalizeData(all.combined.integrated.epithelial.plot)

#Figure
Dotplot_Zhiwei_Version_NoScale(all.combined.integrated.epithelial.plot, genes_to_plot, "plot_annotation") + theme(axis.text.x = element_text(face = "italic", size = 16))
#Fig S7C
carraro <- readRDS("your/path/GSE150674_Seurat_Object.rds")
carraro <- UpdateSeuratObject(carraro)
DimPlot(carraro)

Idents(carraro) <- "type"
carraro <- subset(carraro, idents = "CO") #control cells only
DimPlot(carraro, group.by = "Cell_Type", label = TRUE, pt.size = 0.5)
#Fig S7D
fig7d.goi <- c("TP63" ,"KRT13", "KRT5", "KRT4", "LGR5", "LGR6")
for (i in seq_along(fig7d.goi)) {
  print(FeaturePlot(adult.proxtodistal.rpca.integrated.epi.cleaned, features = fig7d.goi[i], pt.size = 0.5, order = TRUE) & scale_colour_gradientn(colors =c("lightgrey","#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")))
}
#Fig 7E
#Perform Mapping of Adult Proximal Epithelium from Carraro et al., 2021 onto fetal epithelium
DefaultAssay(all.combined.integrated.epithelium) <- "integrated"
DefaultAssay(carraro) <- "RNA"
carraro <- JoinLayers(carraro)
carraro <- NormalizeData(carraro)
common.features <- intersect(rownames(all.combined.integrated.epithelium), rownames(carraro))
length(x = common.features)

all.combined.integrated.epithelium.anchors <- FindTransferAnchors(
  reference = all.combined.integrated.epithelium,
  query = carraro,
  reference.assay = "integrated",
  query.assay = "RNA",
  features = common.features,
  reference.reduction = "pca",
  k.filter = 200, 
)


carraro.mapped.fetal.epithelium <- MapQuery(
  anchorset = all.combined.integrated.epithelium.anchors,
  query = carraro,
  reference = all.combined.integrated.epithelium,
  refdata = "annotation_lvl3",
  reference.reduction = "pca", 
  reduction.model = "umap",
  
)


all.combined.integrated.epithelium$id <- 'reference'
carraro.mapped.fetal.epithelium$id <- 'query'
refquery <- merge(carraro.mapped.fetal.epithelium, all.combined.integrated.epithelium)
refquery[["umap"]] <- merge(carraro.mapped.fetal.epithelium[["ref.umap"]], all.combined.integrated.epithelium[["umap"]])
refquery$predicted.id <- factor(refquery$predicted.id, levels = identities)


DimPlot(refquery, group.by = 'predicted.id', shuffle = TRUE, label = TRUE, cols = color_map) + NoAxes() + NoLegend()
