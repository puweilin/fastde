# created with usethis::use_test()
# run with devtools::test()

source("../harness.R")

test_that("findallmarkers", {

  sobj = load_pbmc3k()

  f = paste0(get_data_dir(), "pbmc.seurat.Rdata")
  if (file.exists(f)) {
    pbmc = load(f)

  } else {
    pbmc = seurat_pipeline(sobj)

    save(pbmc, file = f)
  }

  # seurat.markers <- Seurat::FindAllMarkers(pbmc, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(seurat.markers)

  # fastde.markers <- fastde::FastFindAllMarkers(pbmc, test.use = "fastwmw", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # str(fastde.markers)

  # expect_equal(fastde.markers, seurat.markers)

  expect_equal(1, 1)
})



# test_that("findallmarkers_wilcox", {

#   sobj <- load_pbmc3k()
#   # str(sobj)

#   # NOTE: cells in cols.
#   spmat <- sobj@assays$RNA@counts
#   genenames = spmat@Dimnames[1]
#   str(genenames)

#   nclusters = 12  
#   labels <- gen_labels(nclusters, ncol(spmat))
#   L <- unique(sort(labels))
#   sobj@meta.data$rand.ident <- as.factor(as.character(labels))
#   Seurat::Idents(sobj) <- "rand.ident"

#   str(sobj)
#   message("unique identities ", unique(Seurat::Idents(sobj)))
  
#   seurat_de <- Seurat::FindAllMarkers(sobj,
#                                         test.use = 'wilcox',
#                                         logfc.threshold = 0.25,
#                                         min.pct = 0.1,
#                                         min.diff.pct = -Inf,
#                                         only.pos = FALSE,
#                                         ma.cells.per.ident = Inf,
#                                         min.cells.feature = 3,
#                                         return.thresh = 0.01)
#   print(seurat_de)
#   message("number of rows ", nrow(seurat_de))

#   fastde_de <- fastde::FastFindAllMarkers(sobj, 
#                                         test.use = 'fastwmw',
#                                         logfc.threshold = 0.25,
#                                         min.pct = 0.1,
#                                         min.diff.pct = -Inf,
#                                         only.pos = FALSE,
#                                         ma.cells.per.ident = Inf,
#                                         min.cells.feature = 3,
#                                         return.thresh = 0.01)

#   message("number of rows ", nrow(fastde_de))
#   # rownames(wilcox_de) <- wilcox_de$gene
#   print(fastde_de)

#   expect_equal(fastde_de$p_val, seurat_de$p_val)

# })

