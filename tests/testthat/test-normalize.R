# created with usethis::use_test()
# run with devtools::test()

test_that("normalize_rc", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  # can only handle positive numbers.
  spmat@x = abs(spmat@x)


  seurat_norm <- Seurat::NormalizeData(spmat, normalization.method = "RC", scale.factor=1e4, margin=1, verbose=FALSE)  # margin==1 means across features, i.e. colsum.

  fastde_norm <- fastde::sp_normalize(spmat, normalization.method = "RC", scale.factor=1e4, margin=1, threads=1L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm@x, seurat_norm@x)  # may have small diff due to conversion

  fastde_norm4 <- fastde::sp_normalize(spmat, normalization.method = "RC", scale.factor=1e4, margin=1, threads=4L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm4@x, seurat_norm@x)  # may have small diff due to conversion
})

test_that("normalize_ln", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  # can only handle positive numbers.
  spmat@x = abs(spmat@x)


  seurat_norm <- Seurat::NormalizeData(spmat, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, verbose=FALSE)  # margin==1 means across features, i.e. colsum.

  fastde_norm <- fastde::sp_normalize(spmat, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, threads=1L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm@x, seurat_norm@x)  # may have small diff due to conversion

  fastde_norm4 <- fastde::sp_normalize(spmat, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, threads=4L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm4@x, seurat_norm@x)  # may have small diff due to conversion
})

test_that("normalize_clr", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  # can only handle positive numbers.
  spmat@x = abs(spmat@x)


  seurat_out <- Seurat::NormalizeData(spmat, normalization.method = "CLR", scale.factor=1e4, margin=1, verbose=FALSE)  # margin==1 means across features, i.e. colsum.
  seurat_norm <- as.sparse(x= seurat_out)

  fastde_norm <- fastde::sp_normalize(spmat, normalization.method = "CLR", scale.factor=1e4, margin=1, threads=1L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm@x, seurat_norm@x)  # may have small diff due to conversion

  fastde_norm4 <- fastde::sp_normalize(spmat, normalization.method = "CLR", scale.factor=1e4, margin=1, threads=4L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm4@x, seurat_norm@x)  # may have small diff due to conversion



  seurat_out <- Seurat::NormalizeData(spmat, normalization.method = "CLR", scale.factor=1e4, margin=2, verbose=FALSE)  # margin==1 means across features, i.e. colsum.
  seurat_norm <- as.sparse(x= seurat_out)
  
  fastde_norm <- fastde::sp_normalize(spmat, normalization.method = "CLR", scale.factor=1e4, margin=2, threads=1L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm@x, seurat_norm@x)  # may have small diff due to conversion

  fastde_norm4 <- fastde::sp_normalize(spmat, normalization.method = "CLR", scale.factor=1e4, margin=2, threads=4L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm4@x, seurat_norm@x)  # may have small diff due to conversion
})


test_that("normalize_rc_seurat", {

  sobj <- load_pbmc3k()


  out <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, verbose=FALSE)  # margin==1 means across features, i.e. colsum.
  seurat_norm <- out@assays[[out@active.assay]]@data

  sobj2 <- load_pbmc3k() 
  spmat <- sobj2@assays[[sobj2@active.assay]]@counts

  fastde_norm <- fastde::sp_normalize(spmat, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, threads=1L)  # margin==1 means across features, i.e. colsum.
  expect_equal(fastde_norm@x, seurat_norm@x)  # may have small diff due to conversion

  fastde_norm4 <- fastde::sp_normalize(spmat, normalization.method = "LogNormalize", scale.factor=1e4, margin=1, threads=4L)  # margin==1 means across features, i.e. colsum.

  expect_equal(fastde_norm4@x, seurat_norm@x)  # may have small diff due to conversion
})

