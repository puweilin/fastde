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
#' @import future
seurat_pipeline <- function(sobj) {
    future::plan(sequential)
    
    pbmc <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
    
    
    # vst throws some kind of error or warning: 
    # In simpleLoess(y, x, w, span, degree = degree, parametric = parametric,  :
    # Chernobyl! trL<k 1.8409
    # thi means not enough data poitns to for loess smooth estimation
    pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(pbmc)

    pbmc2 <- Seurat::ScaleData(pbmc, features = all.genes, verbose = FALSE, scale.max = 100000)

    # PCA is using OMP_NUM_THREADS and is not paralleliezd in Seurat
    # if component == 50, irlba fails with   "BLAS/LAPACK routine 'DLASCL' gave error code -4"
    # 30 components is fine.
    pbmc <- Seurat::RunPCA(pbmc2, features = VariableFeatures(object = pbmc), verbose = FALSE, npcs=30)

    pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)

    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)

    # in this pipeline, calling RunUMAP causes failure,
    # due to non finite values in matrix.
    # pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)

    return(pbmc)
}