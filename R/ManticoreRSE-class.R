##' ManticoreRSE
##'
##' Manticore RangedSummarizedExperiment
##'
##' @author Per Unneberg
##'
##' @import SummarizedExperiment
##'
setClass("ManticoreRSE", contains = c("RangedSummarizedExperiment"))


.valid.ManticoreRSE <- function(object) {
    if (length(object)) {
        ## if (!(all(unlist(lapply(assays(object), inherits, "ManticoreDF"))))) {
        ##     txt <- sprintf("All ManticoreRSE assay objects must inherit from ManticoreDF")
        ##     return(txt)
        ## }
    }
    NULL
}


setValidity("ManticoreRSE", .valid.ManticoreRSE)


##' ManticoreRSE
##'
##' ManticoreRSE initialization function
##'
##' @param ... options to be passed to RangedSummarizedExperiment
##' @param window.size Window size
##'
##' @return ManticoreRSE object
##'
##' @author Per Unneberg
##'
ManticoreRSE <- function(..., window.size = integer()) {
    rse <- SummarizedExperiment(...)
    if (!is(rse, "RangedSummarizedExperiment"))
        stop("ManticoreRSE must be setup as a RangedSummarizedExperiment")
    mrse <- as(rse, "ManticoreRSE")
    mrse
}
