library(fastde)
library(tictoc)
library(rhdf5)

source(file = 'utils.R')

# read input
input <- rhdf5::h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/block0_values")
labels_all <- as.vector(rhdf5::h5read("/home/tpan/data/gnw2000/gnw2000_truth.h5", "array/block0_values"))
genenames <- rhdf5::h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/axis1")
samplenames <- rhdf5::h5read("/home/tpan/data/gnw2000/gnw2000.h5", "array/axis0")
# wilcox <- rhdf5::h5read("/home/tpan/build/wave/test-wilcox.h5", "array/block0_values")


colnames(input) <- genenames
rownames(input) <- samplenames

# input <- input[, 1:1000]


# count number of labels
#num_labels = nlevels(labels)

# cat(sprintf("test size:  r %d X c %d.\n", nrow(wilcox), ncol(wilcox)))
cat(sprintf("input size:  r %d X c %d\n", nrow(input), ncol(input)))

labels <- as.integer(1 - labels_all[1:nrow(input)])  # inputs are 0 an 1.   wilcox treat class 0 as other in a formula.   fastde treat 0 as first class.   We need to flip this...
cat(sprintf("Labels: %d \n", length(labels)))

L <- unique(sort(labels))

cat(sprintf("Labels unique: %d \n", length(L)))

# typeof(labels)
# for ( i in labels) {
#     cat(sprintf("%d, ", i))
# }
# cat(sprintf("\n"))


# # print(fastdewilcox_df)
# # fastdewilcox_df$p_val
# # fastdewilcox_df$cluster
# # fastdewilcox_df$genes
# comparemat("c++ vs fastde tstat2", wilcox, fastdewilcox)


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde")
fastdewilcox <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde u")
# # fastdewilcox_stat <- fastde::wmw_fast(input, labels, rtype=as.integer(3), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
# toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde z")
# # fastdewilcox_z <- fastde::wmw_fast(input, labels, rtype=as.integer(4), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
# toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df")
fastdewilcox_df <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

# print(fastdewilcox2_df)


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde4")
fastdewilcox4 <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde u4")
# # fastdewilcox_stat4 <- fastde::wmw_fast(input, labels, rtype=as.integer(3), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
# toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde z4")
# # fastdewilcox_z4 <- fastde::wmw_fast(input, labels, rtype=as.integer(4), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
# toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df4")
fastdewilcox_df4 <- fastde::wmw_fast(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()



cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastdev")
fastdewilcoxv <- fastde::wmw_fastv(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde uv")
# # fastdewilcox_statv <- fastde::wmw_fastv(input, labels, rtype=as.integer(3), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
# toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde zv")
# # fastdewilcox_zv <- fastde::wmw_fastv(input, labels, rtype=as.integer(4), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(1))
# toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_dfv")
fastdewilcox_dfv <- fastde::wmw_fastv(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(1))
toc()

# print(fastdewilcox2_df)


cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde4v")
fastdewilcox4v <- fastde::wmw_fastv(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde u4v")
# # fastdewilcox_stat4v <- fastde::wmw_fastv(input, labels, rtype=as.integer(3), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
# toc()

# cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
# tic("fastde z4v")
# # fastdewilcox_z4v <- fastde::wmw_fastv(input, labels, rtype=as.integer(4), 
#     continuity_correction=TRUE, as_dataframe = FALSE, threads = as.integer(4))
# toc()

cat(sprintf("input %d X %d\n", nrow(input), ncol(input)))
tic("fastde_df4v")
fastdewilcox_df4v <- fastde::wmw_fastv(input, labels, rtype=as.integer(2), 
    continuity_correction=TRUE, as_dataframe = TRUE, threads = as.integer(4))
toc()


# tic("Limma")
# # time and run LIMMA
# Limmawilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
# for ( gene in 1:ncol(input) ) {
#     x <- matrix(as.vector(input[, gene]), ncol=1)
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
# warnings()



# time and run wilcox.test
tic("R builtin")
Rwilcox <- matrix(, ncol = ncol(input), nrow = length(L) )
Rwilcox_stat <- matrix(, ncol = ncol(input), nrow = length(L) )
Rwilcox_z <- matrix(, ncol = ncol(input), nrow = length(L) )
for ( gene in 1:ncol(input) ) {
    x <- matrix(as.vector(input[, gene]), ncol=1)
    # cat(sprintf("gene %d \n", gene))
    # cat(sprintf("x size:  r %d X c %d.\n", nrow(x), ncol(x)))
    i <- 1
    for ( c in L ) {
        # lab <- labels %in% c

        xx <- x[which(labels == c)]
        yy <- x[which(labels != c)]

        # v <- wilcox.test(x ~ lab, alternative="two.sided", correct=TRUE)
        v <- wilcox.test(x = xx, y= yy, alternative="two.sided", correct=TRUE)
        # cat(sprintf("R wilcox %f\n", v))
        Rwilcox[i, gene] <- v$p.value
        # cat(sprintf("R wilcox %f\n", v))
        # Rwilcox_stat[i, gene] <- v$statistic

        # nx = length(xx)
        # ny = length(yy)
        # z <- v$statistic - nx * ny * 0.5
        # r <- rank(c(xx, yy))
        # NTIES <- table(r)
        # NTIES2 <- NTIES[which(NTIES > 1)]
        # SIGMA <- sqrt((nx * ny / 12) *
        #                   ((nx + ny + 1)
        #                    - sum(NTIES^3 - NTIES)
        #                    / ((nx + ny) * (nx + ny - 1))))
        #     CORRECTION <- 0 # sign(z) * 0.5
	    # z <- (z - CORRECTION) / SIGMA
        # Rwilcox_z[i, gene] <- z

        i <- i + 1
    }
}
toc()


# fastdewilcox[, 1]
# Limmawilcox[, 1]
# Rwilcox[, 1]

# fastdewilcox[1, ]
# Limmawilcox[1, ]
# Rwilcox[1, ]

# head(Rwilcox[, 11:20])
# head(fastdewilcox[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
comparemat("R vs fastde wilcox", Rwilcox, fastdewilcox)
comparemat("R vs fastde wilcox4", Rwilcox, fastdewilcox4)
comparemat("R vs fastde wilcoxv", Rwilcox, fastdewilcoxv)
comparemat("R vs fastde wilcox4v", Rwilcox, fastdewilcox4v)
# comparemat("R vs Limma wilcox", Rwilcox, Limmawilcox)


# head(Rwilcox_stat[, 11:20])
# head(fastdewilcox_stat[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
# comparemat("R vs fastde wilcox stat", Rwilcox_stat, fastdewilcox_stat)
# comparemat("R vs fastde wilcox stat4", Rwilcox_stat, fastdewilcox_stat4)
# comparemat("R vs fastde wilcox statv", Rwilcox_stat, fastdewilcox_statv)
# comparemat("R vs fastde wilcox stat4v", Rwilcox_stat, fastdewilcox_stat4v)

# head(Rwilcox_z[, 11:20])
# head(fastdewilcox_z[ , 11:20])

# comparemat("c++ vs R wilcox", wilcox, Rwilcox)
# comparemat("c++ vs fastde wilcox", wilcox, fastdewilcox)
# comparemat("R vs fastde wilcox z", Rwilcox_z, fastdewilcox_z)
# comparemat("R vs fastde wilcox z4", Rwilcox_z, fastdewilcox_z4)
# comparemat("R vs fastde wilcox zv", Rwilcox_z, fastdewilcox_zv)
# comparemat("R vs fastde wilcox z4v", Rwilcox_z, fastdewilcox_z4v)

# pos <- which((Rwilcox_z[1, ] - fastdewilcox_z[1, ]) != 0)
# Rwilcox_z[1, pos[1:10]]
# fastdewilcox_z[1, pos[1:10]]
# # time and run wilcox.test
# tic("R builtin 2sided")
# Rwilcox2 <- matrix(, ncol = ncol(input), nrow = length(L) )
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


# ## compare by calculating the residuals.
# res = Rwilcox - fastdewilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("R naive vs fastde residual tail = %f\n", residual))



# res = Limmawilcox - Rwilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("Limma vs R naive residual = %f\n", residual))

# res = Limmawilcox - fastdewilcox
# residual = sqrt(mean(res * res))

# cat(sprintf("Limma vs fastde residual = %f\n", residual))

