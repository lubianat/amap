.noGenerics <- TRUE
.conflicts.OK <- TRUE

.onLoad <- .First.lib <- function(lib, pkg)
{
    library.dynam("amap", pkg, lib)
    vers <- R.Version()$major
    if(vers <2)
      {
        have.mva <- "package:mva" %in% search()
        if(!have.mva) require("mva")
      }
    else
      {
        have.stats <- "package:stats" %in% search()
        if(!have.stats)   require("stats")
      }

    if("Biobase" %in% .packages(all=TRUE))
      require("Biobase")
    else
        exprs <<- function(x)
            stop("You need to install Biobase package to use this object")
  }

