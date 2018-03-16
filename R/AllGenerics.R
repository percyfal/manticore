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

##' Convert object to GRanges
##'
##' @export
setGeneric("asGRanges", function(x, ...) standardGeneric("asGRanges"))


##' Generic plot function
##'
##' @export
setGeneric("gplot",
           function(data, x=NULL, y=NULL,
                    type="point",
                    xlim=NULL, ylim=NULL, main=NULL,
                    xlab=NULL, ylab=NULL, size=integer(0),
                    colour=list(), colour.var=character(0),
                    wrap=FALSE, wrap.formula=character(0),
                    wrap.ncol=integer(0), compact.facet=TRUE,
                    strip.position="right", scales="free_y",
                    hide.legend=TRUE, hide.xaxis=TRUE, grid=FALSE,
                    text.size=integer(0), text.x.angle=integer(0),
                    text.x.hjust=integer(0), which=NULL, ...)
    standardGeneric("gplot"),
    signature="data")

##' Generic boxplot function
##'
##' @export
setGeneric("gboxplot",
           function(data, formula=character(0),
                    type=character(0),
                    xlim=NULL, ylim=NULL, main=NULL,
                    xlab=NULL, ylab=NULL,
                    colour=list(), colour.var=character(0),
                    wrap=FALSE, wrap.formula=chacarter(0),
                    wrap.ncol=integer(0), compact.facet=TRUE,
                    strip.position="right", scales="free_y",
                    hide.legend=TRUE, hide.xaxis=TRUE, grid=FALSE,
                    text.size=integer(0), text.x.angle=integer(0),
                    text.x.hjust=integer(0), which=NULL, ...)
    standardGeneric("gboxplot"),
    signature="data")
