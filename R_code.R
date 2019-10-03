### R script utilizing Seurat package (v2.0.1) to generate tSNE plots and clustering

library(Seurat)
library(tidyverse)

##-- 10x
#- Load the data
raw_data <- read.csv('raw.expMatrix.csv', header = TRUE, row.names = 1)
dim(raw_data) # 62,167 cells and 20,223 genes - Already filtered

#- Create Seurat object
seuratObj <- CreateSeuratObject(raw.data = raw_data, min.cells = 0, min.genes = 0, project = '10x')
seuratObj # 20,223 genes and 62,167 cells

#- Add meta.data
seuratObj@meta.data$SampleID <- sapply(rownames(seuratObj@meta.data), function(x) str_split(x, '[.]')[[1]][2])
seuratObj@meta.data$Patient <- plyr::mapvalues(x = seuratObj@meta.data$SampleID,
                                               from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                               to = c('CLL-1', 'CLL-1', 'CLL-1', 'CLL-1', 'NHL-9', 'NHL-9', 'NHL-9', 'CLL-2', 'CLL-2', 'CLL-2', 'CLL-2', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-10', 'NHL-9'))
seuratObj@meta.data$Patient <- factor(seuratObj@meta.data$Patient, levels = c('CLL-1', 'CLL-2', 'NHL-9', 'NHL-10'))
seuratObj@meta.data$TimePoint <- plyr::mapvalues(x = seuratObj@meta.data$SampleID,
                                                 from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                                 to = c('IP', 'd21', 'd38', 'd112', 'd12', 'd29', 'd102', 'IP', 'd12', 'd29', 'd83', 'IP', 'd12', 'd28', 'd89', 'IP'))
seuratObj@meta.data$TimePoint <- factor(seuratObj@meta.data$TimePoint, levels = c('IP', 'd12', 'd21', 'd28', 'd29', 'd38', 'd83', 'd89', 'd102', 'd112'))
seuratObj@meta.data$Group <- plyr::mapvalues(x = seuratObj@meta.data$SampleID,
                                             from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                             to = c('IP', 'ExpansionPeak', 'Contraction', 'Late', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP', 'ExpansionPeak', 'Contraction', 'Late', 'IP'))
seuratObj@meta.data$Group <- factor(seuratObj@meta.data$Group, levels = c('IP', 'ExpansionPeak', 'Contraction', 'Late'))
seuratObj@meta.data$Disease <- plyr::mapvalues(x = seuratObj@meta.data$SampleID,
                                               from = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'),
                                               to = c('CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'CLL', 'CLL', 'CLL', 'CLL', 'NHL', 'NHL', 'NHL', 'NHL', 'NHL'))
mito.genes <- grep(pattern = '^MT-', x = rownames(x = seuratObj@data), value = TRUE) # 13 mitochondrial genes
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ]) / Matrix::colSums(seuratObj@raw.data)
seuratObj@meta.data$percent.mito <- percent.mito

#- Normalization
seuratObj <- NormalizeData(object = seuratObj, normalization.method = 'LogNormalize', scale.factor = 10000) # normalized matrix is used as input for MAST analysis

#- Highly variable genes
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
seuratObj@var.genes <- head(rownames(seuratObj@hvg.info)[!str_detect(string = rownames(seuratObj@hvg.info), pattern = '^TR')], 1000) # TRV-genes are removed

#- Dealing with confounders
# nUMI is used as proxy for CDR. Patient is also corrected to avoid any patient effect in clustering and tSNE visualization
seuratObj <- ScaleData(object = seuratObj, do.par = TRUE, num.cores = 6, model.use = 'linear', vars.to.regress = c('nUMI', 'percent.mito', 'Patient'))

#- PCA
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, pcs.compute = 100)
PCElbowPlot(object = seuratObj, num.pc = 100)

#- tSNE
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:15, do.fast = TRUE)

#- Clustering
seuratObj <- FindClusters(object = seuratObj, dims.use = 1:15)

#- Visualization
data_ggplot <- data.table(seuratObj@meta.data,
                          tSNE_1 = seuratObj@dr$tsne@cell.embeddings[, 1],
                          tSNE_2 = seuratObj@dr$tsne@cell.embeddings[, 2])
data_ggplot$Group <- factor(x = data_ggplot$Group, levels = c('IP', 'ExpansionPeak', 'Contraction', 'Late'))
ggplot(data_ggplot) + theme(panel.background = element_rect(colour = "black", size = 1), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), 
                            axis.text.y = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks = element_blank(),
                            legend.text = element_text(size = 10)) + geom_point(aes(x = tSNE_1, y = tSNE_2, color = Group), size = 1, alpha = 1) + labs(x = 'tSNE 1', y = 'tSNE 2', title = '') + guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)), shape = guide_legend(override.aes = list(size = 2, alpha = 1))) + scale_color_discrete('') + scale_shape_discrete('')
data_ggplot <- melt(data_ggplot, id = c('tSNE_1', 'tSNE_2', 'Patient', 'Group'))
ggplot(data_ggplot) + theme(panel.background = element_rect(colour = "black", size = 1), 
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), 
                            axis.text.y = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks = element_blank()) + geom_point(aes(x = tSNE_1, y = tSNE_2, color = Group), size = 1) + labs(x = 'tSNE 1', y = 'tSNE 2') + facet_wrap(~ Patient, ncol = 3) + scale_color_discrete('')
