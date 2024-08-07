#' Quality Control and Normalization of a Seurat Object
#' 
#' This assumes that the input file is already a seurat object in an .rds format.
#' Other files/extensions that is not a seurat object in rds format should be converted first.
#' A package "write_seurat_object" is available on git
#' This package contains three functions: for the quality control \
#' of the seurat object, for the normalization of the seurat object \
#' after QC, and for the visualization of the seurat object during \
#' quality control
#' 
#' @param input is the seurat object with .rds extension
#' @examples 
#' quality_control_seurat("seurat_object.rds") # will output a file after quality filtering the output file should be used as the input for normalization
#' normalization_seurat("afterQC_seurat.rds")
#' @export

qc_seurat <- function(input) {
   seurat_obj <- readRDS(input)
   seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
   seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
   saveRDS(seurat_obj, file = "afterQC_seurat.rds")
}

normalization_seurat <- function(input) {
   seurat_obj <- readRDS(input)
   seurat_obj <- NormalizeData(seurat_obj, normalization.method = "logNormalize", scale.factor = 10000)
   saveRDS(seurat_obj, file = "afternormalization.rds")
}

# this should take in as input the seurat object before QC and normalization
visual_qcmetrices <- function(input){
   seurat_obj <- readRDS(input)
   seurat_obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
   
vplot <- Seurat::VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend() %>%
   FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
   
   ggplot2::ggsave(file = "violin_plot.png")
   
   # to visualzie separate features
   vplot1 <- Seurat::VlnPlot(seurat_obj, features = "nFeature_RNA") + NoLegend()
   vplot2 <- Seurat::VlnPlot(seurat_obj, features = "nCount_RNA") + NoLegend()
   vplot1 <- Seurat::VlnPlot(seurat_obj, features = "percent.mt") + NoLegend()
   
   ggvplot1 <- ggplotly(vplot1)
   ggvplot2 <- ggplotly(vplot2)
   ggvplot3 <- ggplotly(vplot3)
   
annotations = list(list(x = 0.15, y = 1, text = "nFeature_RNA", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE), list(x = 0.5, y = 1, text = "nCount_RNA", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE), list(x = 0.85, y = 1, text = "percent.mt", xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", showarrow = FALSE))

three_plots <- subplot(ggvplot1, ggvplot2, ggvplot3) %>%
   layout(title = "Violin Plots", annotations = annotations)

ggplot2::ggsave(three_plots, file = "violin_plot.png")
}
