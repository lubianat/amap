## Hierarchical clustering, on raw input data; we will use Euclidean
## distance.  A range of criteria are supported; also there is a
## storage-economic option.
##
## We use the general routine, `hc', which caters for 7 criteria,
## using a half dissimilarity matrix; (BTW, this uses the very efficient
## nearest neighbor chain algorithm, which makes this algorithm of
## O(n^2) computational time, and differentiates it from the less
## efficient -- i.e. O(n^3) -- implementations in all commercial
## statistical packages -- as far as I am aware -- except Clustan.)
##
## Clustering Methods:
##
## 1. Ward's minimum variance or error sum of squares method.
## 2. single linkage or nearest neighbor method.
## 3. complete linkage or diameter.
## 4. average linkage, group average, or UPGMA method.
## 5. McQuitty's or WPGMA method.
## 6. median, Gower's or WPGMC method.
## 7. centroid or UPGMC method (7).
##
## Original author: F. Murtagh, May 1992
## R Modifications: Ross Ihaka, Dec 1996
##		    Friedrich Leisch, Apr 1998, Jun 2000
##                  Antoine Lucas, Oct 2002

hclust <-
  function (d, method = "complete", members = NULL) 
  { 
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                   "median", "centroid")
      method <- pmatch(method, METHODS)
      if (is.na(method)) 
        stop("invalid clustering method")
      if (method == -1) 
        stop("ambiguous clustering method")
      n <- as.integer(attr(d, "Size"))
      if (is.null(n)) 
        stop("invalid dissimilarities")
      if (n < 2) 
        stop("Must have n >= 2 objects to cluster")
      labels <- attr(d, "Labels")
      len <- n * (n - 1)/2
      if (is.null(members)) 
        members <- rep(1, n)
      if (length(members) != n) 
        stop("Invalid length of members")
      hcl <- .C("hclust",
                n = as.integer(n),
                len = as.integer(len), 
                method = as.integer(method),
                ia = integer(n),
                ib = integer(n),
                order = integer(n),
                crit = double(n),
                members = as.double(members),
                diss = as.double(d),
                PACKAGE ="amap" )
      
      tree <- list(merge = cbind(hcl$ia[1:(n - 1)],
                     hcl$ib[1:(n - 1)]),
                   height = hcl$crit[1:(n - 1)],
                   order = hcl$order, 
                   labels = attr(d, "Labels"),
                   method = METHODS[method], 
                   call = match.call())
      if (!is.null(attr(d, "method"))) {
        tree$dist.method <- attr(d, "method")
      }
      class(tree) <- "hclust"
      tree
    }


  
plot.hclust <-
    function (x, labels = NULL, hang = 0.1,
              axes = TRUE, frame.plot = FALSE, ann = TRUE,
              main = "Cluster Dendrogram",
              sub = NULL, xlab = NULL, ylab = "Height", ...)
{
    merge <- x$merge
    if (!is.matrix(merge) || ncol(merge) != 2)
	stop("invalid dendrogram")
    n <- nrow(merge)
    height <- as.double(x$height)
    order <- as.double(order(x$order))

    labels <-
	if(missing(labels) || is.null(labels)) {
	    if (is.null(x$labels))
		paste(1:(n+1))
	    else
		as.character(x$labels)
	} else {
	    if(is.logical(labels) && !labels)# FALSE
		character(n+1)
	    else
		as.character(labels)
       }

    plot.new()
    .Internal(dend.window(n, merge, height, order, hang, labels, ...))
    .Internal(dend       (n, merge, height, order, hang, labels, ...))
    if(axes)
        axis(2, at=pretty(range(height)))
    if (frame.plot)
        box(...)
    if (ann) {
        if(!is.null(cl <- x$call) && is.null(sub))
            sub <- paste(deparse(cl[[1]])," (*, \"", x$method,"\")",sep="")
        if(is.null(xlab) && !is.null(cl))
            xlab <- deparse(cl[[2]])
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
    }
    invisible()
}

## For S ``compatibility'': was in cluster just as
## plclust <- plot.hclust ## .Alias

plclust <- function(tree, hang = 0.1, unit = FALSE, level = FALSE, hmin = 0,
                    square = TRUE, labels = NULL, plot. = TRUE,
                    axes = TRUE, frame.plot = FALSE, ann = TRUE,
                    main = "", sub = NULL, xlab = NULL, ylab = "Height")
{
    if(!missing(unit) && unit)		.NotYetUsed("unit", error = FALSE)
    if(!missing(level) && level)	.NotYetUsed("level", error = FALSE)
    if(!missing(hmin) && hmin != 0)	.NotYetUsed("hmin",  error = FALSE)
    if(!missing(square) && !square)	.NotYetUsed("square",error = FALSE)
    if(!missing(plot.) && !plot.)	.NotYetUsed("plot.", error = TRUE)
    plot.hclust(x = tree, labels = labels, hang = hang,
                axes = axes, frame.plot = frame.plot, ann = ann,
                main = main, sub = sub, xlab = xlab, ylab = ylab)
}


