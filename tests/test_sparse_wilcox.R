library(rhdf5)
library(fastde)
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

# read input
# tic("read")
# labels <- as.vector(h5read("/home/tpan/build/wave/labels.h5", "array/block0_values"))
# genenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis1")
# samplenames <- h5read("/home/tpan/build/wave/input.h5", "array/axis0")
# wilcox <- h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")

# input <- rsparsematrix(length(samplenames), length(genenames), 0.1)

# colnames(input) <- genenames
# rownames(input) <- samplenames
# toc()

tic("generate")
# max is 2B.
ncols = 1000  # features  - test with a small number. this is totally parallel
nrows = 10000  # samples
nclusters = 30

clusters = 1:nclusters
labels <- as.integer(sample(clusters, nrows, replace = TRUE))

input <- rsparsematrix(nrows, ncols, 0.05)

samplenames <- as.character(1:nrows);
genenames <- as.character(1:ncols);
colnames(input) <- genenames
rownames(input) <- samplenames
toc()


tic("get unique labels")
# count number of labels
#num_labels = nlevels(labels)

# cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))
cat(sprintf("Labels: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

toc()

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))


tic("Densify")
dinput = as.matrix(input)
toc()

# time and run wilcox.test
tic("R builtin")
Rwilcox <- matrix(, ncol = ncols, nrow = length(L) )
for ( gene in 1:ncol(dinput) ) {
    # x <- matrix(as.vector(input[, gene]), ncol=1)
    x <- dinput[, gene]
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        # lab <- labels %in% c

        v <- wilcox.test(x = x[which(labels == c)], y= x[which(labels != c)], alternative="two.sided", correct=TRUE)
        # cat(sprintf("R wilcox %f\n", v))
        Rwilcox[i, gene] <- v$p.value
        i <- i + 1
    }
}
toc()

# tic("Limma")
# # time and run LIMMA
# Limmawilcox <- matrix(, ncol = ncols, nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     # x <- matrix(as.vector(input[, gene]), ncol=1)
#     x <- dinput[, gene]
#     # cat(sprintf("gene %d \n", gene))
#     # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    
#     i <- 1
#     for ( c in L ) {
#         ind <- which(labels %in% c)
        
#         ## two sided
#         v <- min(2 * min(limma::rankSumTestWithCorrelation(index = ind, statistics = x)), 1)  # two sides.
        
#         ## less
#         # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[1]  # less, left tail
#         ## greater
#         # v <- limma::rankSumTestWithCorrelation(index = ind, statistics = x)[2]    # greater, right tail
#         # cat(sprintf("limma %f\n", v))
#         Limmawilcox[i, gene] <- v
#         i <- i + 1
#     }
# }
# toc()
# # warnings()



cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde")

fastdewilcox <- fastde::wmw_fast(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde_df")
fastdewilcox_df <- fastde::wmw_fast(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()



cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde 4")
fastdewilcox4 <- fastde::wmw_fast(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde_df 4")
fastdewilcox_df4 <- fastde::wmw_fast(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()


cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastdev")
fastdewilcoxv <- fastde::wmw_fastv(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde_dfv")
fastdewilcox_dfv <- fastde::wmw_fastv(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()



cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde 4v")
fastdewilcox4v <- fastde::wmw_fastv(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input %d X %d\n", nrow(dinput), ncol(dinput)))
tic("fastde_df 4v")
fastdewilcox_df4v <- fastde::wmw_fastv(dinput, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()



# print(fastdewilcox_df)
# fastdewilcox_df$p_val
# fastdewilcox_df$cluster
# fastdewilcox_df$genes








# fastdewilcox[, 1]
# Limmawilcox[, 1]
# Rwilcox[, 1]

# fastdewilcox[1, ]
# Limmawilcox[1, ]
# Rwilcox[1, ]


# # time and run wilcox.test
# tic("R builtin 2sided")
# Rwilcox2 <- matrix(, ncol = ncols, nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
#     # cat(sprintf("gene %d \n", gene))
#     # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
#     i <- 1
#     for ( c in L ) {
#         lab <- labels %in% c

#         v <- wilcox.test(x ~ lab, alternative="two.sided")$p.value
#         # cat(sprintf("R wilcox %f\n", v))
#         Rwilcox2[i, gene] <- v
#         i <- i + 1
#     }
# }
# toc()

# Rwilcox2[, 1]


# compare by calculating the residuals.
comparemat("R vs fastde", Rwilcox, fastdewilcox)
# comparemat("limma vs fastde", Limmawilcox, fastdewilcox)

comparemat("R vs fastde", Rwilcox, fastdewilcoxv)
comparemat("R vs fastde", Rwilcox, fastdewilcox4)
comparemat("R vs fastde", Rwilcox, fastdewilcox4v)



#==================



cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde sparse")
fastdewilcox <- fastde::sparse_wmw_fast(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df sparse")
fastdewilcox_df <- fastde::sparse_wmw_fast(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde sparse 4")
fastdewilcox4 <- fastde::sparse_wmw_fast(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df sparse 4")
fastdewilcox_df4 <- fastde::sparse_wmw_fast(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde sparsev")
fastdewilcoxv <- fastde::sparse_wmw_fastv(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df sparsev")
fastdewilcox_dfv <- fastde::sparse_wmw_fastv(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde sparse 4v")
fastdewilcox4v <- fastde::sparse_wmw_fastv(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df sparse 4v")
fastdewilcox_df4v <- fastde::sparse_wmw_fastv(input, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()




## compare by calculating the residuals.
comparemat("R vs sparse fastde", Rwilcox, fastdewilcox)
# comparemat("limma vs sparse fastde", Limmawilcox, fastdewilcox)

comparemat("R vs sparse fastde", Rwilcox, fastdewilcox4)
comparemat("R vs sparse fastde", Rwilcox, fastdewilcoxv)
comparemat("R vs sparse fastde", Rwilcox, fastdewilcox4v)


print("CONVERT To 64bit")
input64 <- as.dgCMatrix64(input)

cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde sparse64")
fastdewilcox64 <- fastde::sparse_wmw_fast(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde_df sparse64")
fastdewilcox_df64 <- fastde::sparse_wmw_fast(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde sparse644")
fastdewilcox644 <- fastde::sparse_wmw_fast(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde_df sparse644")
fastdewilcox_df644 <- fastde::sparse_wmw_fast(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()

cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde sparse64v")
fastdewilcox64v <- fastde::sparse_wmw_fastv(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()


cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde_df sparse64v")
fastdewilcox_df64v <- fastde::sparse_wmw_fastv(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde sparse644v")
fastdewilcox644v <- fastde::sparse_wmw_fastv(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()


cat(sprintf("input64 %d X %d\n", nrow(input64), ncol(input64)))
tic("fastde_df sparse644v")
fastdewilcox_df644v <- fastde::sparse_wmw_fastv(input64, labels, features_as_rows = FALSE, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()



## compare by calculating the residuals.
comparemat("R vs sparse64 fastde", Rwilcox, fastdewilcox64)
# comparemat("limma vs sparse64 fastde", Limmawilcox, fastdewilcox)

comparemat("R vs sparse64 fastde", Rwilcox, fastdewilcox64v)
comparemat("R vs sparse64 fastde", Rwilcox, fastdewilcox644)
comparemat("R vs sparse64 fastde", Rwilcox, fastdewilcox644v)



