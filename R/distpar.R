distpar <- function(x, method="euclidean", nbproc = 2, diag=FALSE, upper=FALSE)
{
    ## account for possible spellings of euclid?an
    if(!is.na(pmatch(method, "euclidian")))
	method <- "euclidean"

    METHODS <- c("euclidean", "maximum",
                 "manhattan", "canberra", "binary","pearson","correlation")
    method <- pmatch(method, METHODS)
    if(is.na(method))
	stop("invalid distance method")
    if(method == -1)
	stop("ambiguous distance method")

    N <- nrow(x <- as.matrix(x))
    d <- .C("R_distancepar",
	    x = as.double(x),
	    nr= N,
	    nc= ncol(x),
	    d = double(N*(N - 1)/2),
	    diag  = as.integer(FALSE),
	    method= as.integer(method),
            nbproc = as.integer(nbproc),
            err = as.integer(0),
	    DUP = FALSE, NAOK=TRUE, PACKAGE="amap")$d
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
