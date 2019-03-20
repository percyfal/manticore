##' summary.ManticoreRSE
##'
##' Summarize the statistics of a ManticoreRSE object
##'
##'
##' @param obj ManticoreRSE object
##' @param which which assay to summarize
##' @param fun summary function
##' @param formula formula to calculate
##' @param per.region normalize score by window size
##' @param ... additional arguments
##'
##' @return a data frame
##' @author Per Unneberg
##' @export
##'
summary.ManticoreRSE <- function(obj, which = NULL, fun = "mean",
                                 per.site = FALSE,
                                 per.region = FALSE,
                                 variable = "score", sites = "sites",
                                 ...) {
    fun <- match.arg(fun, c("mean", "sd", "median", "var"))
    if (is.null(which))
        which <- names(assays(obj))
    else
        which <- match.arg(which, names(assays(obj)), several.ok = TRUE)
    cData <- colData(obj)
    .getData <- function(y, measure) {
        data <- subset(assay(y, measure), select = rownames(cData)[cData$variable == variable])
        if (per.site) {
            sites <- subset(assay(y, measure), select = rownames(cData)[cData$variable == sites])
            sites[sites == 0] <- NA
            data <- data / sites
        }
        data
    }
    y <- do.call("rbind", lapply(which, function(x) {
                              if (!per.region)
                                  apply(.getData(obj, x), 2, fun, ...)
                              else
                                  do.call("rbind", tapply(obj, seqnames(obj), function(y) {apply(.getData(y, x), 2, fun, ...)}))
                          }))
    if (!per.region)
        rownames(y) <- which
    y
}


##' Aggregate ManticoreRSE statistics
##'
##' @param x ManticoreRSE object
##' @param formula formula
##' @param FUN function to apply
##' @param per.site normalize window scores by window length (kb)
##' @param ... additional arguments passed to aggregate
##'
##' @export
##'
setMethod("aggregate", "ManticoreRSE", function(x, formula, FUN, which = NULL, per.site=FALSE, ...) {
    NULL
})
