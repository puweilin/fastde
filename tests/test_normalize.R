library(Seurat)
library(rhdf5)
library(fastde)   # loading library can take a while (3 s) if not preloaded here.
library(tictoc)
library(Matrix)

comparemat <- function(name, A, B) {
    diff <- A - B
    maxdiff <- max(diff)
    mindiff <- min(diff)
    mediandiff <- median(diff)
    meandiff <- mean(diff * diff)
    stdevdiff <- sd(diff * diff)
    cat(sprintf("%s : diff range [%.17g, %.17g], median %.17g, mean %.17g, sd %.17g\n", 
        name, mindiff, maxdiff, mediandiff, meandiff, stdevdiff))
    # mxpos = which(diff == maxdiff, arr.ind = TRUE)
    # mnpos = which(diff == mindiff, arr.ind = TRUE)

    # if ( abs(maxdiff) > .Machine$double.eps)
    #    cat(sprintf("%s : max diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
    #         name, mxpos, A[mxpos], B[mxpos], diff[mxpos]))
    # if  ( abs(mindiff) > .Machine$double.eps)
    #     cat(sprintf("%s : min diff at pos %d:  A %.17g - B %.17g = DIFF %.17g.\n", 
    #         name, mnpos, A[mnpos], B[mnpos], diff[mnpos]))
    
}


# cat(sprintf("NOTE: rounding is done differently in R and C++ (and python):\nR and Python: round half to even. c++: round half away from zero.\n\te.g. R: round(%f) = %f.  C would return %f\n\tMATRIX COMPARISON RESULTS ARE EXPECTED TO BE DIFFERENT BUT WITH ZERO MEAN.\n\tFASTDE FoldChange results are not rounded.\n", -0.5, round(-0.5), -1.0))


tic("generate input")


# rows = 15052
# cols = 38000
rows = 15052
cols = 2800

M <- rsparsematrix(rows, cols, 0.04)
rownames(M) <- paste0("r", 1:rows)
colnames(M) <- paste0("c", 1:cols)

M@x = abs(M@x)

toc()

str(M)


# time and run wilcox.test
tic("Seurat relcount")

seurat_out <- Seurat::NormalizeData(M, normalization.method = "RC", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.

toc()
str(seurat_out)

tic("Fastde relcount")

fastde_out <- fastde::sp_normalize(M, normalization.method = "RC", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.

toc()

str(fastde_out)

comparemat("R vs fastde relcnt ", seurat_out@x, fastde_out@x)



# time and run wilcox.test
tic("Seurat lognorm")

seurat_out <- Seurat::NormalizeData(M, normalization.method = "LogNormalize", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.

toc()
str(seurat_out)

tic("Fastde lognorm")

fastde_out <- fastde::sp_normalize(M, normalization.method = "LogNormalize", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.

toc()

str(fastde_out)

comparemat("R vs fastde ", seurat_out@x, fastde_out@x)



# time and run wilcox.test
tic("Seurat clr")

seurat_out <- Seurat::NormalizeData(M, normalization.method = "CLR", scale.factor=1e4, margin=2)  # margin==1 means across features, i.e. colsum.
seurat_out <- as.sparse(x= seurat_out)

toc()
str(seurat_out)

tic("Fastde clr")

fastde_out <- fastde::sp_normalize(M, normalization.method = "CLR", scale.factor=1e4, margin=2)  # margin==1 means across features, i.e. colsum.

toc()

str(fastde_out)

comparemat("R vs fastde clr ", seurat_out@x, fastde_out@x)



# time and run wilcox.test
tic("Seurat clr row")

seurat_out <- Seurat::NormalizeData(M, normalization.method = "CLR", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.
seurat_out <- as.sparse(x= seurat_out)

toc()
str(seurat_out)

tic("Fastde clr row")

fastde_out <- fastde::sp_normalize(M, normalization.method = "CLR", scale.factor=1e4, margin=1)  # margin==1 means across features, i.e. colsum.

toc()

str(fastde_out)

comparemat("R vs fastde clr row", seurat_out@x, fastde_out@x)
