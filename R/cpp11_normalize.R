
#' R Sparse transpose
#'
#' This implementation directly constructs the new transposed matrix.  
#'     There is random memory writes.
#' 
#' @rdname sp_normalize
#' @param spmat a sparse matrix, of the form dgCMatrix or dgCMatrix64
#' @param normalization.method Method for normalization.
#'  \itemize{
#'   \item{LogNormalize: }{Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. This is then
#'   natural-log transformed using log1p.}
#'   \item{CLR: }{Applies a centered log ratio transformation}
#'   \item{RC: }{Relative counts. Feature counts for each cell are divided by the total
#'   counts for that cell and multiplied by the scale.factor. No log-transformation is applied.
#'   For counts per million (CPM) set \code{scale.factor = 1e6}}
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization
#' @param margin If performing CLR normalization, normalize across features (1) or cells (2)
#' @param threads Number of threads for parallelization
#' @return normalized sparse matrix.
#' @name sp_normalize
#' @concept preprocessing
#' @export
sp_normalize <- function(spmat, 
    normalization.method = 'LogNormalize', 
    scale.factor = 1e4, 
    margin = 1, 
    threads = 1) {
    met <- switch(
        EXPR = normalization.method,
        'LogNormalize' = 0,
        'CLR' = 1,
        'RC' = 2,
        stop("Unknown normalization method: ", normalization.method)
    )
    if (is(spmat, 'dgCMatrix')) {
        
        newx <- cpp11_sp_normalize(x=spmat@x, i=spmat@i, p=spmat@p, nrow=spmat@Dim[1], ncol=spmat@Dim[2], scale_factor=scale.factor, margin=margin, method=met, threads=threads)
        # str(newx)
        out <- new("dgCMatrix", x=newx, i=spmat@i, p=spmat@p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        return (out)
    } else if (is(spmat, 'dgCMatrix64')) {
        newx <- cpp11_sp64_normalize(x=spmat@x, i=spmat@i, p=spmat@p, nrow=spmat@Dim[1], ncol=spmat@Dim[2], scale_factor=scale.factor, margin=margin, method=met, threads=threads)
        # str(newx)
        out <- new("dgCMatrix64", x=newx, i=spmat@i, p=spmat@p, Dim=spmat@Dim, Dimnames=list(rownames(spmat), colnames(spmat)))
        return(out)
    } else {
        print("ERROR: unsupported data type for normalize")
        return(spmat)
    }
}

