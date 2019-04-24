##' GENOMEList
##'
##' @param obj an R object
##' @param ... additional arguments
##'
##' @rdname GENOMEList
##' @export
##'
setGeneric("GENOMEList", function(obj, ...) standardGeneric("GENOMEList"))

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
