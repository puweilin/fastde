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
sp_normalize <- function(spmat, normalization.method = 'LogNormalize', scale.factor = 1e4, margin = 1, threads = 1) {
    met <- switch(
        EXPR = normalization.method,
        'LogNormalize' = 0,
        'CLR' = 1,
        'RC' = 2,
        stop("Unknown normalization method: ", normalization.method)
    )
    if (is(spmat, 'dgCMatrix')) {
        
        tic("normalize")
        newx <- cpp11_sp_normalize(x=spmat@x, i=spmat@i, p=spmat@p, nrow=spmat@Dim[1], ncol=spmat@Dim[2], scale_factor=scale.factor, margin=margin, method=met, threads=threads)
        toc()
        str(newx)
        tic("create obj")
        out <- new("dgCMatrix", x=newx, i=spmat@i, p=spmat@p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        toc()
        return (out)
    } else if (is(spmat, 'dgCMatrix64')) {
        tic("normalize")
        newx <- cpp11_sp64_normalize(x=spmat@x, i=spmat@i, p=spmat@p, nrow=spmat@Dim[1], ncol=spmat@Dim[2], scale_factor=scale.factor, margin=margin, method=met, threads=threads)
        toc()
        str(newx)
        tic("create obj")
        out <- new("dgCMatrix64", x=newx, i=spmat@i, p=spmat@p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        toc()
        return(out)
    } else {
        print("USING R DEFAULT")
        return(spmat)
    }
}

