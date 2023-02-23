library(tictoc)

#' R Sparse transpose
#'
#' This implementation directly constructs the new transposed matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_transpose
#' @param spmat a sparse matrix, of the form dgCMatrix
#' @return matrix dense matrix.
#' @name sp_transpose
#' @export
sp_normalize <- function(spmat, method = 'LogNormalize', scale_factor = 0.5, threads = 1) {
    met <- switch(
        EXPR = method,
        'LogNormalize' = 0,
        'CLR' = 1,
        'RC' = 2,
        stop("Unknown normalization method: ", method)
    )
    if (is(spmat, 'dgCMatrix')) {
        
        tic("normalize")
        mlist <- cpp11_sp_normalize(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], met, scale_factor, threads)
        toc()
        tic("create obj")
        out <- new("dgCMatrix", x=mlist$x, i=mlist$i, p=mlist$p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        toc()
        return (out)
    } else if (is(spmat, 'dgCMatrix64')) {
        tic("normalize")
        mlist <- cpp11_sp64_normalize(spmat@x, spmat@i, spmat@p, spmat@Dim[1], spmat@Dim[2], met, scale_factor, threads)
        toc()
        tic("create obj")
        out <- new("dgCMatrix64", x=mlist$x, i=mlist$i, p=mlist$p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        toc()
        return(out)
    } else {
        print("USING R DEFAULT")
        return(t(spmat))
    }
}

