.noGenerics <- TRUE
.conflicts.OK <- TRUE

.onLoad <- .First.lib <- function(lib, pkg)
{
    library.dynam("amap", pkg, lib)
    have.mva <- "package:mva" %in% search()
    if(!have.mva) require("mva")
}

