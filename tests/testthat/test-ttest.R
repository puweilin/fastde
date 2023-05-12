# created with usethis::use_test()
# run with devtools::test()

test_that("dense_ttest", {

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

  Rttest <- matrix(, ncol = ncols, nrow = length(L) )
  for ( gene in 1:ncols ) {
      dat <- as.vector(input[, gene])
      # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
      # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
      i <- 1
      for ( c in L ) {
          lab <- labels %in% c

          v <- t.test(x = dat[lab], y = dat[!lab])$p.value
          # cat(sprintf("R wilcox %f\n", v))
          Rttest[i, gene] <- v
          i <- i + 1
      }
  }
  rownames(Rttest) <- as.character(L)
  colnames(Rttest) <- colnames(input)

  fastdettest <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, 
    threads = 1L, alternative = as.integer(2), var_equal = FALSE)

  expect_equal(Rttest, fastdettest)  # may have small diff due to conversion

  fastdettest4 <- fastde::ttest_fast(input, labels, as_dataframe = FALSE, 
    threads = 4L, alternative = as.integer(2), var_equal = FALSE)

  expect_equal(Rttest, fastdettest4)  # may have small diff due to conversion
})


test_that("sparse_ttest", {

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
  Rttest <- matrix(, ncol = ncols, nrow = length(L) )
  for ( gene in 1:ncol(input) ) {
      dat <- as.vector(input[, gene])
      # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
      # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
      i <- 1
      for ( c in L ) {
          lab <- labels %in% c

          v <- t.test(x = dat[lab], y = dat[!lab])$p.value
          # cat(sprintf("R wilcox %f\n", v))
          Rttest[i, gene] <- v
          i <- i + 1
      }
  }
  rownames(Rttest) <- as.character(L)
  colnames(Rttest) <- colnames(input)

  fastdettest <- fastde::sparse_ttest_fast(spmat, labels,
    features_as_rows = FALSE,  as_dataframe = FALSE, threads = 1L, alternative =  as.integer(2), var_equal = FALSE)


  expect_equal(Rttest, fastdettest)  # may have small diff due to conversion

  fastdettest4 <- fastde::sparse_ttest_fast(spmat, labels,
    features_as_rows = FALSE,  as_dataframe = FALSE, threads = 4L, alternative =  as.integer(2), var_equal = FALSE)

  expect_equal(Rttest, fastdettest4)  # may have small diff due to conversion
})

test_that("sparse_ttest_64", {

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
  Rttest <- matrix(, ncol = ncols, nrow = length(L) )
  for ( gene in 1:ncol(input) ) {
      dat <- as.vector(input[, gene])
      # if ( (gene %% 100) == 0 ) cat(sprintf("R building t.test.  gene %d \n", gene))
      # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
      i <- 1
      for ( c in L ) {
          lab <- labels %in% c

          v <- t.test(x = dat[lab], y = dat[!lab])$p.value
          # cat(sprintf("R wilcox %f\n", v))
          Rttest[i, gene] <- v
          i <- i + 1
      }
  }
  rownames(Rttest) <- as.character(L)
  colnames(Rttest) <- colnames(input)

  spmat64 = as.dgCMatrix64(spmat)

  fastdettest <- fastde::sparse_ttest_fast(spmat64, labels,
    features_as_rows = FALSE,  as_dataframe = FALSE, threads = 1L, alternative =  as.integer(2), var_equal = FALSE)


  expect_equal(Rttest, fastdettest)  # may have small diff due to conversion

  fastdettest4 <- fastde::sparse_ttest_fast(spmat64, labels,
    features_as_rows = FALSE,  as_dataframe = FALSE, threads = 4L, alternative =  as.integer(2), var_equal = FALSE)

  expect_equal(Rttest, fastdettest4)  # may have small diff due to conversion
})

