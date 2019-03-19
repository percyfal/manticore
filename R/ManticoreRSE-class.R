.valid.ManticoreRSE <- function(object) {
    message("Validating ManticoreRSE")
    message(lapply(assays(object), inherits, "ManticoreDF"))
    if (length(object)) {
        if (!(all(unlist(lapply(assays(object), inherits, "ManticoreDF")))))
            message("All ManticoreRSE assay objects must inherit from ManticoreDF")
    }
}

##' ManticoreRSE
##'
##' Manticore RangedSummarizedExperiment
##'
##' @param object
##' @return
##' @author Per Unneberg
##'
##' @import SummarizedExperiment
##'
setClass("ManticoreRSE",
         contains = c("RangedSummarizedExperiment"))

setValidity2("ManticoreRSE", .valid.ManticoreRSE)


ManticoreRSE <- function(...) {
    rse <- SummarizedExperiment(...)
    if (!is(rse, "RangedSummarizedExperiment"))
        stop("ManticoreRSE must be setup as a RangedSummarizedExperiment")
    as(rse, "ManticoreRSE")
}
