# quality-control-seurat
This package contains two functions a) `qc_seurat()` and b) `normalization_seurat()`.

The first function performs quality control on a Seurat object generated from scRNAseq files. This removes cells with a low number of unique genes that imply dying or low-quality cells, cells or likely doublets/multiplets with high gene count, and cells with a high percentage of mitochondrial genes which imply they are dying cells.

The second function is performed after quality control. It normalizes the feature expression measurements for each cell by the total expression then multiplied by a scale factor and log-transforms the result.

## Installation
The package can be installed using
```
devtools::install_github("BTIP/quality-control-seurat")
```

## Example
The output of this function would be an rds file which when loaded in R would be a Seurat object.
```
# to perform qc
qc_seurat("seurat_object.rds")

# to normalize the data
normalization_seurat("afterQC_seurat.rds")
```

## scRNAseq processing workflow
The standard scRNAseq processing workflow with the R package Seurat consists of seven (7) steps. The output of this package and function should be used as input for the scRNAseq processing pipeline. 

The following are the repositories of the packages for every step of the pipeline:
1. QC and filtering: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
2. Normalization: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
3. Identification of highly variable features: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
4. Scaling: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
5. Linear Dimensionality Reduction (PCA): [pcaseurat package](https://github.com/BTIP2024/pca-seurat)
6. Clustering: [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)
7. Non-linear dimensionality reduction (t-SNE and UMAP): [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)

An overview of the pipeline and its outputs can be observed below:
![](https://github.com/user-attachments/assets/de8e812a-6c6d-475a-9b58-33156348de11)
