##' @rdname window.size
##' @description Get window.size from a WindowedSummarizedExperiment object
##'
##' @param obj WindowedSummarizedExperiment object
##'
setMethod("window.size", "WindowedSummarizedExperiment",
          function(obj) {
    return (rowRanges(obj)@window.size)
})
