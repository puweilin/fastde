# TODO: [] limit row and col dimensions to int.     need for long per axis?
# NOTE: be extremely careful about changing the tags for documentation.

## annotations for documentation can be tricky.  particularly around
## signature of generic methods.  see 
## https://stackoverflow.com/questions/7356120/how-to-properly-document-s4-methods-using-roxygen2

# do not include @name in the distached methos
# https://github.com/r-lib/roxygen2/issues/695


#' Based on dgCMatrix class (definitions.R) to support compressed sparse column matrice with
#' up to 2^53 non-zero entries.   R's dgCMatrix allows 2^31 or 2 billion non-zero elements.
#' 
#' @slot i integer array of row ids.
#' @slot p floating array of positions offsets.  support 2^53 non-zeros based on IEEE double's support for exact integers
#' @importClassesFrom Matrix dsparseMatrix
# ' @name dgCMatrix64
# ' @rdname dgCMatrix64
# ' @aliases dgCMatrix64-class
# ' @exportClass dgCMatrix64
dgCMatrix64 <- setClass("dgCMatrix64", contains = c("dsparseMatrix"),
	 slots = c(i = "integer", p = "numeric"),
	 prototype = prototype(i = 0L, p = 0.0)
     )
# no validation.  use Dim - assume dimensions are at most 2B..

# use initialize
setMethod("initialize", "dgCMatrix64",
    function(.Object, x,
                     i,
                     p,
                     Dim,
                     Dimnames=c(NULL, NULL)
                     )  
{
    .Object <- callNextMethod()

    .Object@x <- x
    .Object@i <- as.integer(i)
    .Object@p <- as.double(p)
    .Object@Dim <- as.integer(Dim)
    .Object@Dimnames <- Dimnames

    .Object
}
)

# # create new dgCMatrix64
# .newdgCMatrix64 <-  function(x=NULL,
#                      i=NULL,
#                      p=0,
#                      Dim=c(0, 0),
#                      Dimnames=c(NULL, NULL)
#                      )
# {
#     message("newdgCMatrix64 being called...")
#     newx <- new("dgCMatrix64")
#     slot(newx,"x",check=FALSE) <- x
#     # x <- NULL
#     slot(newx,"i",check=FALSE) <- as.integer(i)
#     # i <- NULL
#     slot(newx,"p",check=FALSE) <- as.double(p)
#     # p <- NULL
#     slot(newx,"Dim",check=FALSE) <- as.integer(Dim)
#     # Dim <- NULL
#     slot(newx,"Dimnames", check=FALSE) <- Dimnames
#     # Dimnames <- NULL
    
#     return(newx)
# }

#' Check if an object is a dgCMatrix64 object

#' @param obj object to test
#' @return True if obj is a dgCMatrix64 object
#' @docType methods
#' @rdname is.dgCMatrix64-methods
#' @name is.dgCMatrix64
#' @export
is.dgCMatrix64 <- function(obj) 
    is(obj, "dgCMatrix64")

# #' @exportMethod = as.dgCMatrix64
#' @name as.dgCMatrix64
#' @docType methods
#' @param obj an object
#' @param eps value type conversion tolerance, unused
#' @return a dgCMatrix64 object 
#' @rdname as.dgCMatrix64-methods
#' @export
as.dgCMatrix64 <- function(obj, eps = .Machine$double.eps)
    stop("coercion not defined form this class")

# # ' convert a dgCMatrix64 object to another dgCMatrix64 object.

#' @docType methods
#' @rdname as.dgCMatrix64-methods
#' @param obj a dgcMatrix64 object
#' @param eps value type conversion tolerance, unused
#' @return a dgCMatrix64 object
#' @aliases as.dgCMatrix64,dgCMatrix64,ANY-method
#' @export
#' @method as.dgCMatrix64 dgCMatrix64
as.dgCMatrix64.dgCMatrix64 <- function(obj, eps = .Machine$double.eps)  {
    newobj <- new("dgCMatrix64", x = obj@x, i = obj@i, p = obj@p,
        Dim = obj@Dim, Dimnames = obj@Dimnames)
    newobj
}

#' convert a dgCMatrix object to another dgCMatrix64 object.
#' 
#' @docType methods
#' @param obj a dgcMatrix64 object
#' @param eps value type conversion tolerance, unused
#' @return a dgCMatrix64 object 
#' @method as.dgCMatrix64 dgCMatrix
#' @rdname as.dgCMatrix64-methods
#' @aliases as.dgCMatrix64,dgCMatrix,ANY-method
#' @export
as.dgCMatrix64.dgCMatrix <- function(obj, eps = .Machine$double.eps)  {
    newobj <- new("dgCMatrix64", x = obj@x, i = obj@i, p = obj@p,
        Dim = obj@Dim, Dimnames = obj@Dimnames)
    newobj
}

# do this in c++ (sparsify)
##' @export 
#as.dgCMatrix64.matrix <- function(obj, eps = .Machine$double.eps) {

#' @name as.dgCMatrix64
#' @param obj object to be converted to dgCMatrix64
#' @param eps error tolerance for data type conversion.  unused.
#' @return dgCMatrix64 sparse matrix
#' @docType methods
#' @rdname as.dgCMatrix64-methods
#' @export
setGeneric("as.dgCMatrix64", as.dgCMatrix64)
# setGeneric("as.dgCMatrix64", function(obj, eps = .Machine$double.eps)  {
#     standardGeneric("as.dgCMatrix64")
# })


#' @param obj object to be converted to dgCMatrix64
#' @param eps error tolerance for data type conversion.  unused.
#' @return dgCMatrix64 sparse matrix
#' @docType methods
#' @aliases as.dgCMatrix64,dgCMatrix64,ANY-method
#' @rdname as.dgCMatrix64-methods
#' @export
setMethod("as.dgCMatrix64", "dgCMatrix64", as.dgCMatrix64.dgCMatrix64)
# setMethod("as.dgCMatrix64", signature("dgCMatrix64"), function(obj, eps = .Machine$double.eps)  {
#     new('dgCMatrix64', x = obj@x, i = obj@i, p = obj@p, Dim=obj@Dim, Dimnames = obj@Dimnames)
# })


#' @param obj object to be converted to dgCMatrix64
#' @param eps error tolerance for data type conversion.  unused.
#' @return dgCMatrix64 sparse matrix
#' @docType methods
#' @aliases as.dgCMatrix64,dgCMatrix,ANY-method
#' @rdname as.dgCMatrix64-methods
#' @export
setMethod("as.dgCMatrix64", "dgCMatrix", as.dgCMatrix64.dgCMatrix)
# setMethod("as.dgCMatrix64", signature("dgCMatrix"), function(obj, eps = .Machine$double.eps)  {
#     new('dgCMatrix64', x = obj@x, i = obj@i, p = as.numeric(obj@p), Dim=obj@Dim, Dimnames = obj@Dimnames)
# })

#setMethod("as.dgCMatrix64","matrix", as.dgCMatrix64.matrix)

