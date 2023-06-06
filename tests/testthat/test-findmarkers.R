# created with usethis::use_test()
# run with devtools::test()

source("../harness.R")

test_that("findallmarkers_wilcox", {

  sobj = load_pbmc3k()

  f = paste0(get_data_dir(), "pbmc.seurat.Rdata")
  if (file.exists(f)) {
    load(f)   # this replaces the variables of the same names in the workspace.

  } else {
    # str(sobj)
    pbmc = seurat_pipeline(sobj)

    save(pbmc, file = f)
  }
  # str(pbmc)

  seurat.markers <- Seurat::FindAllMarkers(pbmc, test.use = "wilcox", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(seurat.markers)

  fastde.markers <- fastde::FastFindAllMarkers(pbmc, test.use = "fastwmw", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(fastde.markers)

  expect_equal(fastde.markers$cluster, seurat.markers$cluster)
  expect_equal(fastde.markers$p_val, seurat.markers$p_val)
  expect_equal(fastde.markers$avg_log2FC, seurat.markers$avg_log2FC)
  expect_equal(fastde.markers$p_val_adj, seurat.markers$p_val_adj)
  expect_equal(round(fastde.markers$pct.1, digits = 3), seurat.markers$pct.1)
  expect_equal(round(fastde.markers$pct.2, digits = 3), seurat.markers$pct.2)
  # expect_equal(fastde.markers$gene, seurat.markers$gene)
})


test_that("findallmarkers_ttest", {

  sobj = load_pbmc3k()

  f = paste0(get_data_dir(), "pbmc.seurat.Rdata")
  if (file.exists(f)) {
    load(f)   # this replaces the variables of the same names in the workspace.

  } else {
    # str(sobj)
    pbmc = seurat_pipeline(sobj)

    save(pbmc, file = f)
  }
  # str(pbmc)

  seurat.markers <- Seurat::FindAllMarkers(pbmc, test.use = "t", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(seurat.markers)

  fastde.markers <- fastde::FastFindAllMarkers(pbmc, test.use = "fast_t", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(fastde.markers)

  expect_equal(fastde.markers$cluster, seurat.markers$cluster)
  expect_equal(fastde.markers$p_val, seurat.markers$p_val)
  expect_equal(fastde.markers$avg_log2FC, seurat.markers$avg_log2FC)
  expect_equal(fastde.markers$p_val_adj, seurat.markers$p_val_adj)
  expect_equal(round(fastde.markers$pct.1, digits = 3), seurat.markers$pct.1)
  expect_equal(round(fastde.markers$pct.2, digits = 3), seurat.markers$pct.2)
  # expect_equal(fastde.markers$gene, seurat.markers$gene)
})

