##' Create GStats object
##'
##' Retrieve genome stats for
##'
##' @param object An R object
##' @param gr A GRanges object or missing
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
                    gr=NULL,
                    statistics=character(0),
                    use.population.names=FALSE,
                    use.region.names=FALSE,
                    quiet=FALSE,
                    ...)
    standardGeneric("GStats"))


##' GENOMEList
##'
##' @param obj an R object
##' @param ... additional arguments
##'
##' @rdname GENOMEList
##' @export
##'
setGeneric("GENOMEList", function(obj, ...) standardGeneric("GENOMEList"))

##' Convert object to GRanges
##'
##' @param x an R object
##' @param ... additional arguments
##'
##' @rdname asGRanges
##' @export
##'
setGeneric("asGRanges", function(x, ...) standardGeneric("asGRanges"))

##' sites
##'
##' Get the number of sites for an object
##'
##' @param obj A data object, typically containing coordinates and
##'     number of sites per coordinate
##' @return A Vector of sites
##' @author Per Unneberg
##'
##' @export
##' @rdname sites
##'
setGeneric("sites", function(obj) standardGeneric("sites"))

##' Generic plot function
##'
##' @export
##' @rdname gplot
##'
setGeneric("gplot",
           function(data, x=NULL, y=NULL,
                    type="point",
                    xlim=NULL, ylim=NULL, main=NULL,
                    xlab=NULL, ylab=NULL, size=integer(0),
                    colour=list(), colour.var=character(0),
                    wrap=FALSE, wrap.formula=character(0),
                    pos=character(0),
                    wrap.ncol=integer(0), compact.facet=TRUE,
                    strip.position="right", scales="free_y",
                    hide.legend=TRUE, hide.xaxis=TRUE, grid=FALSE,
                    text.size=integer(0), text.x.angle=integer(0),
                    text.x.hjust=integer(0), which=NULL, per.site=FALSE,
                    zscore=FALSE, ...)
    standardGeneric("gplot"),
    signature="data")

##' Generic boxplot function
##'
##' @export
##' @rdname gboxplot
##'
setGeneric("gboxplot",
           function(data, formula=character(0),
                    type=character(0),
                    xlim=NULL, ylim=NULL, main=NULL,
                    xlab=NULL, ylab=NULL,
                    colour=list(), colour.var=character(0),
                    wrap=FALSE, wrap.formula=chacarter(0),
                    wrap.ncol=integer(0), compact.facet=TRUE,
                    strip.position="right", scales="free_y",
                    hide.legend=TRUE, grid=FALSE,
                    text.size=integer(0), text.x.angle=integer(0),
                    text.x.hjust=integer(0), which=NULL, ...)
    standardGeneric("gboxplot"),
    signature="data")


##' Generic window.size function
##'
##' @export
##' @rdname window.size
##'
setGeneric("window.size",
           function(obj)
    standardGeneric("window.size"),
    signature="obj")
