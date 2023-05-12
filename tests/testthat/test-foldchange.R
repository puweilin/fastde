# created with usethis::use_test()
# run with devtools::test()

# foldchange in Seurat are rounded to 3 decimal places.
# expect_equal uses all.equal, which specified tolerance as a relative vale
# but we need absolute.

# round fastde output instead.

test_that("dense_foldchange", {

  nrows = 3000
  ncols = 10
  nclusters = 12
  spmat <- rsparsematrix(nrows, ncols, 0.05)

  samplenames <- as.character(1:nrows);
  genenames <- as.character(1:ncols);
  colnames(spmat) <- genenames
  rownames(spmat) <- samplenames

  input <- as.matrix(spmat)
  
  labels = gen_labels(nclusters, nrows)
  L <- unique(sort(labels))  # depending on nrows, may be fewer than nclusters
  clusternames = as.character(L)  

  seuratfc <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc1 <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc2 <- matrix(, ncol = ncols, nrow = length(L) )
  
  dimnames(seuratfc) <- list(clusternames, genenames)

  # need features as rows.
  x <- t(input)
  colnames(x) <- samplenames
  rownames(x) <- genenames

  for ( c in L ) {
      cells.1 <- which(labels %in% c)
      cells.2 <- which(! (labels %in% c))

      v <- Seurat::FoldChange(x, cells.1, cells.2, mean.fxn=rowMeans, fc.name="fc")

      # cat(sprintf("R wilcox %f\n", v))
      pos <- which(L == c)
      seuratfc[pos, ] <- as.vector(v$fc)
      seuratperc1[pos, ] <- as.vector(v$pct.1)
      seuratperc2[pos, ] <- as.vector(v$pct.2)

      # seuratfc
  }

  fastdefc <- fastde::ComputeFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 1L)

  expect_equal(fastdefc$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion

  fastdefc4 <- fastde::ComputeFoldChange(input, labels, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 4L)

  expect_equal(fastdefc4$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion

})


test_that("sparse_foldchange", {

  nrows = 3000
  ncols = 10
  nclusters = 12
  spmat <- rsparsematrix(nrows, ncols, 0.05)

  samplenames <- as.character(1:nrows);
  genenames <- as.character(1:ncols);
  colnames(spmat) <- genenames
  rownames(spmat) <- samplenames

  labels = gen_labels(nclusters, nrows)
  L <- unique(sort(labels))
  clusternames = as.character(L)

  input = as.matrix(spmat)
  seuratfc <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc1 <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc2 <- matrix(, ncol = ncols, nrow = length(L) )
  rownames(seuratfc) <- clusternames
  colnames(seuratfc) <- genenames

  # need features as rows.
  x <- t(input)
  colnames(x) <- samplenames
  rownames(x) <- genenames

  for ( c in L ) {
      cells.1 <- which(labels %in% c)
      cells.2 <- which(! (labels %in% c))

      v <- Seurat::FoldChange(x, cells.1, cells.2, mean.fxn=rowMeans, fc.name="fc")

      # cat(sprintf("R wilcox %f\n", v))
      pos <- which(L == c)
      seuratfc[pos, ] <- as.vector(v$fc)
      seuratperc1[pos, ] <- as.vector(v$pct.1)
      seuratperc2[pos, ] <- as.vector(v$pct.2)

      # seuratfc
  }


  fastdefc <- fastde::ComputeFoldChangeSparse(spmat, labels, features_as_rows = FALSE, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 1L)

  expect_equal(fastdefc$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion

  fastdefc4 <- fastde::ComputeFoldChangeSparse(spmat, labels, features_as_rows = FALSE, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 4L)

  expect_equal(fastdefc4$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion
})

test_that("sparse_foldchange_64", {

  nrows = 3000
  ncols = 10
  nclusters = 12
  spmat <- rsparsematrix(nrows, ncols, 0.05)

  samplenames <- as.character(1:nrows);
  genenames <- as.character(1:ncols);
  colnames(spmat) <- genenames
  rownames(spmat) <- samplenames

  labels = gen_labels(nclusters, nrows)
  L <- unique(sort(labels))
  clusternames = as.character(L)

  input = as.matrix(spmat)
  seuratfc <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc1 <- matrix(, ncol = ncols, nrow = length(L) )
  seuratperc2 <- matrix(, ncol = ncols, nrow = length(L) )
  rownames(seuratfc) <- clusternames
  colnames(seuratfc) <- genenames

  # need features as rows.
  x <- t(input)
  colnames(x) <- samplenames
  rownames(x) <- genenames

  for ( c in L ) {
      cells.1 <- which(labels %in% c)
      cells.2 <- which(! (labels %in% c))

      v <- Seurat::FoldChange(x, cells.1, cells.2, mean.fxn=rowMeans, fc.name="fc")

      # cat(sprintf("R wilcox %f\n", v))
      pos <- which(L == c)
      seuratfc[pos, ] <- as.vector(v$fc)
      seuratperc1[pos, ] <- as.vector(v$pct.1)
      seuratperc2[pos, ] <- as.vector(v$pct.2)

      # seuratfc
  }

  spmat64 = as.dgCMatrix64(spmat)

  fastdefc <- fastde::ComputeFoldChangeSparse(spmat64, labels, features_as_rows = FALSE, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 1L)

  expect_equal(fastdefc$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion

  fastdefc4 <- fastde::ComputeFoldChangeSparse(spmat64, labels, features_as_rows = FALSE, calc_percents = TRUE, fc_name = "fc", 
    use_expm1 = FALSE, min_threshold = 0.0, use_log = FALSE, log_base = 2.0, 
    use_pseudocount = FALSE, as_dataframe = FALSE, threads = 4L)

  expect_equal(fastdefc4$fc, seuratfc)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.1, digits = 3), seuratperc1)  # may have small diff due to conversion
  expect_equal(round(fastdefc4$pct.2, digits = 3), seuratperc2)  # may have small diff due to conversion
})


