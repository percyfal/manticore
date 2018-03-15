##' Create GStats object
##'
##' Retrieve genome stats for 
##' 
##' @param object An R object
##' @param statistics statistics to return
##' @param use.population.names use population names of object for plotting labels
##' @param use.region.names Use the region names as row names
##' @param quiet suppress messages
##' @param ... Arguments to pass to data access functions
##' @return A GStats object
##' @author Per Unneberg
##'
##' @export
##' @rdname GStats
##' 
setGeneric("GStats",
           function(object,
                    statistics=character(0),
                    use.population.names=FALSE,
                    use.region.names=FALSE,
                    quiet=FALSE,
                    ...)
    standardGeneric("GStats"))


##' GENOMEList
##'
##' Convert objects to GENOMEList
##' @title GENOMEList
##' @param obj passed object
##' @param ... additional arguments
##' @return GENOMEList instance
##' @author Per Unneberg
##' @export
setGeneric("GENOMEList", function(obj, ...) standardGeneric("GENOMEList"))
