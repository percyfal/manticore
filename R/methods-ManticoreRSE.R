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


## override cbind; for combining objects with different ranges
## See https://gist.github.com/PeteHaitch/8993b096cfa7ccd08c13 for discussion
setMethod("cbind", "ManticoreRSE",
          function(..., deparse.level = 1) {
    args <- unname(list(...))
    .merge <- function(...) {
        merge(..., all = TRUE)
    }
    rr <- do.call(.merge, lapply(args, rowRanges))

    newArgs <- lapply(args, function(x){
        data <- lapply(assays(x), function(y) {
            types <- sapply(y, class)
            n <- length(rr)
            z <- DataFrame(lapply(types, function(x) {call(x, n)}))
            z[1:n,] <- NA
            colnames(z) <- colnames(y)
            i <- match(rowRanges(x), rr)
            z[i[!is.na(i)],] <- y
            z})
        SummarizedExperiment(rowRanges = rr, colData = colData(x), assays = data)
    })
    sexp <- do.call(cbind, newArgs)
    ManticoreRSE(rowRanges = rowRanges(sexp), colData = colData(sexp), assays = assays(sexp), window.size = rr@window.size)
})


## For plotting; select assay and column and convert to data fram together with rowRanges
setMethod("as.data.frame", "ManticoreRSE",
          function(x, row.names = NULL, optional = FALSE,
                   long = FALSE,
                   ...) {
    assayData <- as.data.frame(assays(x))
    colnames(assayData)[1:2] <- c("group", "assay")
    if (long)
        assayData <- gather(assayData, select = c(-group, -assay))
    cbind(rowRanges(x), assayData)
})
