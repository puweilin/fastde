# created with usethis::use_test()
# run with devtools::test()

test_that("dense_wilcox", {

  nrows = 3000
  ncols = 10
  nclusters = 12
  spmat <- rsparsematrix(nrows, ncols, 0.05)

  samplenames <- as.character(1:nrows);
  genenames <- as.character(1:ncols);
  colnames(spmat) <- genenames
  rownames(spmat) <- samplenames

  input <- as.matrix(spmat)
  
  labels = gen_labels(nclusters, nrow(input))
  L <- unique(sort(labels))

  Rwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
  for ( gene in 1:ncol(input) ) {
      x <- matrix(as.vector(input[, gene]), ncol=1)
      i <- 1
      for ( c in L ) {

          xx <- x[which(labels == c)]
          yy <- x[which(labels != c)]

          v <- wilcox.test(x = xx, y= yy, alternative="two.sided", correct=TRUE)
          Rwilcox[i, gene] <- v$p.value

          i <- i + 1
      }
  }
  rownames(Rwilcox) <- as.character(L)
  colnames(Rwilcox) <- colnames(input)

  fastdewilcox <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))

  expect_equal(Rwilcox, fastdewilcox)  # may have small diff due to conversion

  fastdewilcox4 <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))

  expect_equal(Rwilcox, fastdewilcox4)  # may have small diff due to conversion
})


test_that("sparse_wilcox", {

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


  input = as.matrix(spmat)
  Rwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
  for ( gene in 1:ncol(input) ) {
      x <- matrix(as.vector(input[, gene]), ncol=1)
      i <- 1
      for ( c in L ) {

          xx <- x[which(labels == c)]
          yy <- x[which(labels != c)]

          v <- wilcox.test(x = xx, y= yy, alternative="two.sided", correct=TRUE)
          Rwilcox[i, gene] <- v$p.value

          i <- i + 1
      }
  }
  rownames(Rwilcox) <- as.character(L)
  colnames(Rwilcox) <- colnames(input)

  fastdewilcox <- fastde::sparse_wmw_fast(spmat, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1), features_as_rows = FALSE)

  expect_equal(Rwilcox, fastdewilcox)  # may have small diff due to conversion

  fastdewilcox4 <- fastde::sparse_wmw_fast(spmat, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4), features_as_rows = FALSE)

  expect_equal(Rwilcox, fastdewilcox4)  # may have small diff due to conversion
})

test_that("sparse_wilcox_64", {

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


  input = as.matrix(spmat)
  Rwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
  for ( gene in 1:ncol(input) ) {
      x <- matrix(as.vector(input[, gene]), ncol=1)
      i <- 1
      for ( c in L ) {

          xx <- x[which(labels == c)]
          yy <- x[which(labels != c)]

          v <- wilcox.test(x = xx, y= yy, alternative="two.sided", correct=TRUE)
          Rwilcox[i, gene] <- v$p.value

          i <- i + 1
      }
  }
  rownames(Rwilcox) <- as.character(L)
  colnames(Rwilcox) <- colnames(input)

  spmat64 = as.dgCMatrix64(spmat)

  fastdewilcox <- fastde::sparse_wmw_fast(spmat64, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1), features_as_rows = FALSE)

  expect_equal(Rwilcox, fastdewilcox)  # may have small diff due to conversion

  fastdewilcox4 <- fastde::sparse_wmw_fast(spmat64, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4), features_as_rows = FALSE)

  expect_equal(Rwilcox, fastdewilcox4)  # may have small diff due to conversion
})

