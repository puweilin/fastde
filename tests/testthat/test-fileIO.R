# created with usethis::use_test()
# run with devtools::test()

test_that("roundtrip", {
  sobj <- load_pbmc3k()

  spmat <- sobj@assays$RNA@counts

  fastde::Write10X_h5(spmat, "pbmc3k_spmat.h5")

  spmat2 <- Seurat::Read10X_h5("pbmc3k_spmat.h5")
  
  expect_identical(spmat@x, spmat2@x)
  expect_identical(spmat@i, spmat2@i)
  expect_equal(as.numeric(spmat@p), spmat2@p)  # may have small diff due to conversion
  expect_identical(spmat@Dim, spmat2@Dim)
  expect_identical(spmat@Dimnames, spmat2@Dimnames)
})

test_that("roundtrip32_64", {
  sobj <- load_pbmc3k()

  spmat <- sobj@assays$RNA@counts

  fastde::Write10X_h5(spmat, "pbmc3k_spmat.h5")

  spmat2 <- fastde::Read10X_h5_big("pbmc3k_spmat.h5")
  
  expect_identical(spmat@x, spmat2@x)
  expect_identical(spmat@i, spmat2@i)
  expect_equal(as.numeric(spmat@p), spmat2@p)  # may have small diff due to conversion
  expect_identical(spmat@Dim, spmat2@Dim)
  expect_identical(spmat@Dimnames, spmat2@Dimnames)
})


test_that("roundtrip64", {
  sobj <- load_pbmc3k()

  spmat <- sobj@assays$RNA@counts
  spmat64 <- as.dgCMatrix64(spmat)

  fastde::Write10X_h5(spmat64, "pbmc3k_spmat64.h5")

  spmat2 <- fastde::Read10X_h5_big("pbmc3k_spmat64.h5")
  
  expect_identical(spmat@x, spmat2@x)
  expect_identical(spmat@i, spmat2@i)
  expect_equal(as.numeric(spmat@p), spmat2@p)  # may have small diff due to conversion
  expect_identical(spmat@Dim, spmat2@Dim)
  expect_identical(spmat@Dimnames, spmat2@Dimnames)
})

