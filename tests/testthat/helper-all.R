# if (!require("SeuratData")) {
#     devtools::install_github('satijalab/seurat-data')
# }
# library(SeuratData)


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
    datadir = system.file("testdata", package = "fastde")

    return(datadir)
}

# create small pbmc3k dataset
## do not call during test.
make_small_pbmc3k <- function(count) {
    installed <- intersect(x = paste0("pbmc3k", ".", "SeuratData"), y = rownames(x = SeuratData::InstalledData()))
    if (length(x = installed) == 0) {
        SeuratData::InstallData("pbmc3k")
    }
    data("pbmc3k")

    # now sample.  tried to used "variable features" and it resulted in most features having counts of 1 or 2.
    # where original data has 258 levels.
    # so instead, just pick the first n features.  exclude NA in the annotations
    spbmc <- pbmc[1:2000, !is.na(pbmc@meta.data$seurat_annotations)]

    # Don't bother setting the existing seurat annotation as labels.  just rerun clustering.
    # Seurat::Idents(spbmc) <- spbmc@meta.data$seurat_annotations

    #save to disk for later reload.
    fastde::Write10X_h5(spbmc@assays$RNA@counts, paste0(get_data_dir(), "/pbmc3k_2k.h5"))
}


##' @importFrom SeuratData InstallData
#' @import Seurat
load_pbmc3k <- function() {
    # tried 1.) packaging up the pbmc3k data - files too big.
    #       2.) SeuratData - missing Matrix package during check()
    #       3.) TENxPBMCData - need to convert to seurat obj.  This is not available with r-devel...
    #       4.) subsample 1 ahead of time and use the saved data.


    # 3  # using Bioconductor TENxPBMCData as input
    # pbmc3k.sce=TENxPBMCData('pbmc3k')
    # # must convert from delayed matrix to sparse matrix.
    # counts(pbmc3k.sce, withDimnames=FALSE) <- as(counts(pbmc3k.sce, withDimnames = FALSE), "dgCMatrix")
    # # now convert to seurat object.  since there is no logcount, data=NULL
    # pbmc <- as.Seurat(pbmc3k.sce, counts = 'counts', data = NULL)
    # # slots are not consistently named. -  use active.assay

    # 2  not part of CRAN and there is version mismatch.
    # installed <- intersect(x = paste0("pbmc3k", ".", "SeuratData"), y = rownames(x = SeuratData::InstalledData()))
    # if (length(x = installed) == 0) {
    #     SeuratData::InstallData("pbmc3k")
    # }
    # data("pbmc3k")
    # # # str(pbmc3k)
    # return(pbmc3k)

    # 1. datadir = get_data_dir()
    # dataset = "pbmc3k"
    # # Load the PBMC dataset
    # spmat <- Seurat::Read10X(data.dir = paste0(datadir, '/', dataset))
    # # Initialize the Seurat object with the raw (non-normalized data).
    # pbmc <- Seurat::CreateSeuratObject(counts = spmat, project = dataset, min.cells = 3, min.features = 200)

    # 4. 
    datadir = get_data_dir()
    dataset = "pbmc3k"
    # Load the PBMC dataset
    spmat <- Seurat::Read10X_h5(paste0(datadir, '/', dataset, '_2k.h5'))
    # Initialize the Seurat object with the raw (non-normalized data).  50 features because of the subsampling.
    pbmc <- Seurat::CreateSeuratObject(counts = spmat, project = dataset, min.cells = 3, min.features = 50)

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