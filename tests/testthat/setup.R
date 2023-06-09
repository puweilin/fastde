library(tictoc)
library(Seurat)
library(SeuratData)
library(Matrix)
library(future)
library(fastde)

SeuratData::InstallData("pbmc3k")
data("pbmc3k")
#str(pbmc3k)

