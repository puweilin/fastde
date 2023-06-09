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
    # dir_path = getwd()  

    # datadir = paste0(dir_path, '/../data/')

    datadir = system.file("testdata", package = "fastde")

    return(datadir)
}


##' @importFrom SeuratData InstallData
#' @import Seurat
load_pbmc3k <- function() {
    # # using Bioconductor TENxPBMCData as input
    # pbmc3k.sce=TENxPBMCData('pbmc3k')
    # # must convert from delayed matrix to sparse matrix.
    # counts(pbmc3k.sce, withDimnames=FALSE) <- as(counts(pbmc3k.sce, withDimnames = FALSE), "dgCMatrix")
    # # now convert to seurat object.  since there is no logcount, data=NULL
    # pbmc <- as.Seurat(pbmc3k.sce, counts = 'counts', data = NULL)
    # slots are not consistent.

    # not part of CRAN and there is version mismatch.
    # SeuratData::InstallData("pbmc3k")
    # data("pbmc3k")
    # str(pbmc3k)
    # return(pbmc3k)

    # alternative - data in packge.
    datadir = get_data_dir()
    dataset = "pbmc3k"
    # Load the PBMC dataset
    spmat <- Seurat::Read10X(data.dir = paste0(datadir, '/', dataset))
    # Initialize the Seurat object with the raw (non-normalized data).
    pbmc <- Seurat::CreateSeuratObject(counts = spmat, project = dataset, min.cells = 3, min.features = 200)

    return(pbmc)


}

#' @import Seurat
#' @import future
seurat_pipeline <- function(sobj) {
    future::plan(sequential)
    
    pbmc <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose=FALSE)
    
    
    # vst throws some kind of error or warning: 
    # In simpleLoess(y, x, w, span, degree = degree, parametric = parametric,  :
    # Chernobyl! trL<k 1.8409
    # this means not enough data points to for loess smooth estimation
    # use all genes for this test.
    #pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
    all.genes <- rownames(pbmc)

    pbmc <- Seurat::ScaleData(pbmc, features = all.genes, verbose = FALSE, scale.max = 100000)

    # PCA is using OMP_NUM_THREADS and is not paralleliezd in Seurat
    # if component == 50, irlba fails with   "BLAS/LAPACK routine 'DLASCL' gave error code -4"
    # 30 components is fine.
    # pbmc <- Seurat::RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE, npcs=30)
    # continuing - use all genes instead of just variable features.
    pbmc <- Seurat::RunPCA(pbmc, features = all.genes, verbose = FALSE, npcs=30)

    pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:30, verbose=FALSE)

    pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5, verbose=FALSE)

    # str(pbmc)

    # in this pipeline, calling RunUMAP causes failure,
    # due to non finite values in matrix.
    # pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)

    return(pbmc)
}