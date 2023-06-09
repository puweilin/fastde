# created with usethis::use_test()
# run with devtools::test()

test_that("findallmarkers_wilcox", {

  # f = paste0(get_data_dir(), "pbmc.seurat.Rdata")
  # if (file.exists(f)) {
  #   load(f)   # this replaces the variables of the same names in the workspace.

  # } else {
    sobj = load_pbmc3k()
  #   # str(sobj)
    pbmc = seurat_pipeline(sobj)

  #   save(pbmc, file = f)
  # }
  # str(pbmc)

  # pbmc = pbmc3k


  seurat.markers <- Seurat::FindAllMarkers(pbmc, test.use = "wilcox", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  fastde.markers <- fastde::FastFindAllMarkers(pbmc, test.use = "fastwmw", verbose=FALSE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(fastde.markers)

  expect_equal(sort(fastde.markers$cluster), sort(seurat.markers$cluster))
  expect_equal(sort(fastde.markers$p_val), sort(seurat.markers$p_val))
  expect_equal(sort(fastde.markers$avg_log2FC), sort(seurat.markers$avg_log2FC))
  expect_equal(sort(fastde.markers$p_val_adj), sort(seurat.markers$p_val_adj))
  expect_equal(sort(round(fastde.markers$pct.1, digits = 3)), sort(seurat.markers$pct.1))
  expect_equal(sort(round(fastde.markers$pct.2, digits = 3)), sort(seurat.markers$pct.2))
  # expect_equal(fastde.markers$gene, seurat.markers$gene)
})


test_that("findallmarkers_ttest", {

#   f = paste0(get_data_dir(), "pbmc.seurat.Rdata")
#   if (file.exists(f)) {
#     load(f)   # this replaces the variables of the same names in the workspace.

#   } else {
    sobj = load_pbmc3k()
#     # str(sobj)
    pbmc = seurat_pipeline(sobj)

#     save(pbmc, file = f)
#   }
#   str(pbmc)
#   pbmc = pbmc3k

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

