library(tictoc)
library(Seurat)

# comparemat <- function(name, A, B) {
#     diff <- A - B
#     maxdiff <- max(diff)
#     mindiff <- min(diff)
#     mediandiff <- median(diff)
#     meandiff <- mean(diff * diff)
#     stdevdiff <- sd(diff * diff)
#     cat(sprintf("%s : diff range [%.17g, %.17g], median %.17g, mean %.17g, sd %.17g\n", 
#         name, mindiff, maxdiff, mediandiff, meandiff, stdevdiff))

#         # mxpos = which(diff == maxdiff, arr.ind = TRUE)
#     # mnpos = which(diff == mindiff, arr.ind = TRUE)

#     # if ( abs(maxdiff) > .Machine$double.eps)
#     #    cat(sprintf("%s : max diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
#     #         name, mxpos, A[mxpos], B[mxpos], diff[mxpos]))
#     # if  ( abs(mindiff) > .Machine$double.eps)
#     #     cat(sprintf("%s : min diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
#     #         name, mnpos, A[mnpos], B[mnpos], diff[mnpos]))
    
# }


print_mat <- function(mat, count) {
    print(head(mat, n = count)[, 1:count])
    print("...")
    print(tail(mat, n = count)[, 1:count])
}


gen_labels <- function(nclusters, nsamples) {
    
    # generate fake class labels.
    clusters = 1:nclusters
    labels <- as.integer(sample(clusters, nsamples, replace = TRUE))

    return(labels)
}

get_data_dir <- function() {
    # should be in the "tests/testthat" directory
    dir_path = getwd()  

    datadir = paste0(dir_path, '/../data/')

    return(datadir)
}


##' @importFrom SeuratData InstallData
#' @import Seurat
load_pbmc3k <- function() {
    datadir = get_data_dir()
    dataset = "pbmc3k"

    # Load the PBMC dataset
    pbmc.data <- Seurat::Read10X(data.dir = paste0(datadir, dataset))
    # Initialize the Seurat object with the raw (non-normalized data).
    pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = dataset, min.cells = 3, min.features = 200)

    # SeuratData::InstallData("pbmc3k")
    # data("pbmc3k.final")
    # return(pbmc3k.final)

    return(pbmc)
}

#' @import Seurat
seurat_pipeline <- function(sobj) {
    pbmc <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
    # future::plan(sequential)
    
    pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(pbmc)

    pbmc <- Seurat::ScaleData(pbmc, features = all.genes, verbose = FALSE, scale.max = 1000000)


    # PCA is using OMP_NUM_THREADS and is not paralleliezd in Seurat
    pbmc <- Seurat::RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE, npcs=50)

    pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)

    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)

    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)

    return(pbmc)
}