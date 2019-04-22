##' WindowedSummarizedExperiment
##'
##' The WindowedSummarizedExperiment class extends the
##' RangedSummarizedExperiment class and adds a requirement that the
##' rowRanges are windows of specified length. The rowRanges object is
##' a Windows object which is an extension of GRanges.
##'
##' @author Per Unneberg
##'
##' @import SummarizedExperiment
##'
setClass("WindowedSummarizedExperiment", contains = c("RangedSummarizedExperiment"))


.valid.WindowedSummarizedExperiment <- function(object) {
    if (length(object)) {
        if (!inherits(rowRanges(object), "Windows")) {
            txt <- sprintf("The rowRanges of a WindowedSummarizedExperiment must inherit from the Windows class")
            return(txt)
        }
    }
    NULL
}


setValidity("WindowedSummarizedExperiment", .valid.WindowedSummarizedExperiment)


##' Windowedsummarizedexperiment
##'
##' WindowedSummarizedExperiment initialization function
##'
##' @param ... options to be passed to RangedSummarizedExperiment
##'
##' @return WindowedSummarizedExperiment object
##'
##' @author Per Unneberg
##'
WindowedSummarizedExperiment <- function(...) {
    rse <- SummarizedExperiment(...)
    if (!is(rse, "RangedSummarizedExperiment"))
        stop("WindowedSummarizedExperiment must be setup as a RangedSummarizedExperiment")
    wse <- as(rse, "WindowedSummarizedExperiment")
    wse
}
