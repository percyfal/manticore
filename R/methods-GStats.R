##' summary.GStats
##'
##' Summarize the statistics of a GStats object by regions or
##' populations
##'
##' @param gs GStats object
##' @param per.site calculate per site values
##' @param fun summary function
##' @param per.region normalize score by window size
##' @param ... additional arguments
##' @return a data frame
##' @author Per Unneberg
##' @export
##'
summary.GStats <- function(gs, per.site=TRUE, fun="mean", per.region=FALSE, ...) {
    fun <- match.arg(fun, c("mean", "sd", "median", "var"))
    df <- assay(gs)
    if (per.site) df <- df / width(gs)
    if (!per.region) {
        m <- matrix(apply(df, 2, fun, ...), ncol = length(levels(gs$statistic)))
        ret <- as.data.frame(m, row.names = levels(gs$population))
        colnames(ret) <- levels(gs$statistic)
    } else {
        m <- do.call("rbind", by(t(df), gs$statistic, function(x) {apply(x, 2, fun, ...)}))
        ret <- as.data.frame(t(m), row.names = SummarizedExperiment::rowData(gs)$feature_id)
    }
    ret
}


##' Aggregate GStats statistics
##'
##' @param x GStats object
##' @param formula formula
##' @param FUN function to apply
##' @param per.site normalize window scores by window length
##' @param ... additional arguments passed to aggregate
##'
##' @export
##'
setMethod("aggregate", "GStats", function(x, formula, FUN, per.site=FALSE, ...) {
    df <- as.data.frame(asGRanges(x, per.site = per.site))
    res <- aggregate(formula, df, FUN, ...)
    res
})

##' Convert GStats to GRanges
##'
##' Convert GStats object to GRanges with values filled by data assay
##' values
##'
##' @param long return long format
##' @param per.site calculate per site statistics (divide by width)
##' @return GRanges instance
##' @author Per Unneberg
##'
##' @rdname asGRanges
##' @export
##'
setMethod("asGRanges", "GStats", function(x, long=TRUE, per.site=FALSE) {
    gr <- SummarizedExperiment::rowRanges(x)
    y <- as.data.frame(assay(x))
    colnames(y) <- make.names(rownames(SummarizedExperiment::colData(x)))
    if (per.site) y <- y / width(x)
    values(gr) <- cbind(as.data.frame(values(gr)), y)
    if (long) {
        df <- tidyr::gather(as.data.frame(values(gr)), statistic, value, -"feature_id")
        i <- match(df$statistic, make.names(rownames(SummarizedExperiment::colData(x))))
        if (any(is.na(i)))
            stop("names in long data frame do not match those of rownames(colData(x)); check for special characters in population names")
        n.rep <- nrow(df) / length(gr)
        gr <- GRanges(
            ranges = IRanges(start = rep(start(gr), n.rep),
                             end = rep(end(gr), n.rep)),
            seqnames = rep(seqnames(gr), n.rep),
            strand = rep(strand(gr), n.rep),
            seqinfo = seqinfo(gr),
            seqlengths = seqlengths(gr))
        values(gr) <- df
        gr$statistic <- x$statistic[i]
        gr$population <- x$population[i]
    }
    gr
})
