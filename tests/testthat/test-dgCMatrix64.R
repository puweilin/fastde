# created with usethis::use_test()
# run with devtools::test()

test_that("create", {
  sobj <- load_pbmc3k()

  spmat <- sobj@assays[[sobj@active.assay]]@counts

  x <- spmat@x
  i <- spmat@i
  p <- spmat@p
  dim <- spmat@Dim
  dimnames <- spmat@Dimnames
  c <- dim[2]
  nelem <- p[c+1]

  dat64 = new('dgCMatrix64', x = x, i = i, p = p, Dim = dim, Dimnames = dimnames)

  expect_vector(dat64@x, ptype = double(), size = nelem)
  expect_vector(dat64@i, ptype = integer(), size = nelem)
  expect_vector(dat64@p, ptype = double(), size = c + 1)
  
  expect_identical(x, dat64@x)
  expect_identical(i, dat64@i)
  expect_equal(as.numeric(p), dat64@p)  # may have small diff due to conversion
  expect_identical(dim, dat64@Dim)
  expect_identical(dimnames, dat64@Dimnames)
})


test_that("as.dgCMatrix64", {
  sobj <- load_pbmc3k()

  spmat <- sobj@assays[[sobj@active.assay]]@counts

  x <- spmat@x
  i <- spmat@i
  p <- spmat@p
  dim <- spmat@Dim
  dimnames <- spmat@Dimnames
  c <- dim[2]
  nelem <- p[c+1]

  dat64 = as.dgCMatrix64(spmat)

  expect_vector(dat64@x, ptype = double(), size = nelem)
  expect_vector(dat64@i, ptype = integer(), size = nelem)
  expect_vector(dat64@p, ptype = double(), size = c + 1)
  
  expect_identical(x, dat64@x)
  expect_identical(i, dat64@i)
  expect_equal(as.numeric(p), dat64@p)  # may have small diff due to conversion
  expect_identical(dim, dat64@Dim)
  expect_identical(dimnames, dat64@Dimnames)
})
