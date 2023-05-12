# created with usethis::use_test()
# run with devtools::test()

test_that("sp_transpose", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  tmat = t(spmat)

  ftmat = fastde::sp_transpose(spmat, threads = 1L)

  expect_identical(ftmat, tmat)

  ftmat4 = fastde::sp_transpose(spmat, threads = 4L)

  expect_identical(ftmat4, tmat)


  spmat64 <- as.dgCMatrix64(spmat)

  ftmat = fastde::sp_transpose(spmat64, threads = 1L)

  expect_identical(ftmat@x, tmat@x)
  expect_identical(ftmat@i, tmat@i)
  expect_equal(ftmat@p, tmat@p)
  expect_identical(ftmat@Dim, tmat@Dim)
  expect_identical(ftmat@Dimnames, tmat@Dimnames)

  ftmat4 = fastde::sp_transpose(spmat64, threads = 4L)

  expect_identical(ftmat4@x, tmat@x)
  expect_identical(ftmat4@i, tmat@i)
  expect_equal(ftmat4@p, tmat@p)
  expect_identical(ftmat4@Dim, tmat@Dim)
  expect_identical(ftmat4@Dimnames, tmat@Dimnames)

})



test_that("sp_to_dense", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  tmat = as.matrix(spmat)

  ftmat = fastde::sp_to_dense(spmat, threads = 1L)

  expect_identical(ftmat, tmat)

  ftmat4 = fastde::sp_to_dense(spmat, threads = 4L)

  expect_identical(ftmat4, tmat)


  spmat64 <- as.dgCMatrix64(spmat)


  ftmat = fastde::sp_to_dense(spmat, threads = 1L)

  expect_identical(ftmat, tmat)

  ftmat4 = fastde::sp_to_dense(spmat, threads = 4L)

  expect_identical(ftmat4, tmat)

})

test_that("sp_to_dense_transposed", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  tmat = as.matrix(t(spmat))

  ftmat = fastde::sp_to_dense_transposed(spmat, threads = 1L)

  expect_identical(ftmat, tmat)

  ftmat4 = fastde::sp_to_dense_transposed(spmat, threads = 4L)

  expect_identical(ftmat4, tmat)


  spmat64 <- as.dgCMatrix64(spmat)

  ftmat = fastde::sp_to_dense_transposed(spmat64, threads = 1L)

  expect_identical(ftmat, tmat)

  ftmat4 = fastde::sp_to_dense_transposed(spmat64, threads = 4L)

  expect_identical(ftmat4, tmat)

})


test_that("sp_rbind", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  spmats <- list(spmat, spmat, spmat)

  tmat = rbind(spmat, spmat, spmat)

  ftmat = fastde::sp_rbind(spmats, threads = 1L, method = 0L)

  expect_identical(ftmat@x, tmat@x)
  expect_identical(ftmat@i, tmat@i)
  expect_equal(ftmat@p, tmat@p)
  expect_identical(ftmat@Dim, tmat@Dim)
  expect_identical(ftmat@Dimnames, tmat@Dimnames)

  ftmat4 = fastde::sp_rbind(spmats, threads = 4L, method = 0L)

  expect_identical(ftmat4@x, tmat@x)
  expect_identical(ftmat4@i, tmat@i)
  expect_equal(ftmat4@p, tmat@p)
  expect_identical(ftmat4@Dim, tmat@Dim)
  expect_identical(ftmat4@Dimnames, tmat@Dimnames)



  spmat64 <- as.dgCMatrix64(spmat)
  spmats64 <- list(spmat64, spmat64, spmat64)

  ftmat = fastde::sp_rbind(spmats64, threads = 1L, method = 0L)

  expect_identical(ftmat@x, tmat@x)
  expect_identical(ftmat@i, tmat@i)
  expect_equal(ftmat@p, tmat@p)
  expect_identical(ftmat@Dim, tmat@Dim)
  expect_identical(ftmat@Dimnames, tmat@Dimnames)

  ftmat4 = fastde::sp_rbind(spmats64, threads = 4L, method = 0L)

  expect_identical(ftmat4@x, tmat@x)
  expect_identical(ftmat4@i, tmat@i)
  expect_equal(ftmat4@p, tmat@p)
  expect_identical(ftmat4@Dim, tmat@Dim)
  expect_identical(ftmat4@Dimnames, tmat@Dimnames)


})

test_that("sp_cbind", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)

  spmats <- list(spmat, spmat, spmat)

  tmat = cbind(spmat, spmat, spmat)

  ftmat = fastde::sp_cbind(spmats, threads = 1L, method = 0L)

  expect_identical(ftmat@x, tmat@x)
  expect_identical(ftmat@i, tmat@i)
  expect_equal(ftmat@p, tmat@p)
  expect_identical(ftmat@Dim, tmat@Dim)
  expect_identical(ftmat@Dimnames, tmat@Dimnames)

  ftmat4 = fastde::sp_cbind(spmats, threads = 4L, method = 0L)

  expect_identical(ftmat4@x, tmat@x)
  expect_identical(ftmat4@i, tmat@i)
  expect_equal(ftmat4@p, tmat@p)
  expect_identical(ftmat4@Dim, tmat@Dim)
  expect_identical(ftmat4@Dimnames, tmat@Dimnames)



  spmat64 <- as.dgCMatrix64(spmat)
  spmats64 <- list(spmat64, spmat64, spmat64)

  ftmat = fastde::sp_cbind(spmats64, threads = 1L, method = 0L)

  expect_identical(ftmat@x, tmat@x)
  expect_identical(ftmat@i, tmat@i)
  expect_equal(ftmat@p, tmat@p)
  expect_identical(ftmat@Dim, tmat@Dim)
  expect_identical(ftmat@Dimnames, tmat@Dimnames)

  ftmat4 = fastde::sp_cbind(spmats64, threads = 4L, method = 0L)

  expect_identical(ftmat4@x, tmat@x)
  expect_identical(ftmat4@i, tmat@i)
  expect_equal(ftmat4@p, tmat@p)
  expect_identical(ftmat4@Dim, tmat@Dim)
  expect_identical(ftmat4@Dimnames, tmat@Dimnames)


})

test_that("sp_colsum", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)


  tmat = colSums(spmat)

  ftmat = fastde::sp_colSums(spmat, threads = 1L, method = 0L)

  expect_equal(ftmat, tmat)

  ftmat4 = fastde::sp_colSums(spmat, threads = 4L, method = 0L)

  expect_equal(ftmat4, tmat)


  spmat64 <- as.dgCMatrix64(spmat)

  ftmat = fastde::sp_colSums(spmat64, threads = 1L, method = 0L)

  expect_equal(ftmat, tmat)

  ftmat4 = fastde::sp_colSums(spmat64, threads = 4L, method = 0L)

  expect_equal(ftmat4, tmat)

})

test_that("sp_colsum", {

  nrows = 3000
  ncols = 10

  spmat <- rsparsematrix(nrows, ncols, 0.05)
  rownames(spmat) <- paste0("r", 1:nrows)
  colnames(spmat) <- paste0("c", 1:ncols)


  tmat = rowSums(spmat)

  ftmat = fastde::sp_rowSums(spmat, threads = 1L, method = 0L)

  expect_equal(ftmat, tmat)

  ftmat4 = fastde::sp_rowSums(spmat, threads = 4L, method = 0L)

  expect_equal(ftmat4, tmat)


  spmat64 <- as.dgCMatrix64(spmat)

  ftmat = fastde::sp_rowSums(spmat64, threads = 1L, method = 0L)

  expect_equal(ftmat, tmat)

  ftmat4 = fastde::sp_rowSums(spmat64, threads = 4L, method = 0L)

  expect_equal(ftmat4, tmat)

})
