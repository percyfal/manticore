##' summary.GStats
##'
##' Summarize the statistics of a GStats object by regions or
##' populations
##'
##' @param gr GRanges object
##' @param per.site calculate per site values
##' @param ... additional arguments
##' @return a data frame
##' @author Per Unneberg
##' @export
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
        ret <- as.data.frame(t(m), row.names = rowData(gs)$feature_id)
    }
    ret
}

##' Convert GStats to GRanges
##'
##' Convert GStats object to GRanges with values filled by data assay
##' values
##'
##' @param x GStats object
##' @param long return long format
##' @param per.site calculate per site statistics (divide by width)
##' @return GRanges instance
##' @author Per Unneberg
##' @export
##'
setMethod("asGRanges", "GStats", function(x, long=TRUE, per.site=FALSE) {
    gr <- rowRanges(x)
    y <- as.data.frame(assay(x))
    colnames(y) <- make.names(rownames(colData(x)))
    if (per.site) y <- y / width(x)
    values(gr) <- cbind(as.data.frame(values(gr)), y)
    if (long) {
        df <- gather(as.data.frame(values(gr)), statistic, value, -"feature_id")
        i <- match(df$statistic, make.names(rownames(colData(x))))
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

## ##' Aggregate statistics over regions in a GStats object
## ##'
## ##'
## ##'
## ##' @rdname aggregate_region_stats
## ##' @export
## setMethod("aggregate_region_stats", "GStats", function(object, agg.fun=c("sum", "mean"), ...) {
##     agg.res <- do.call("rbind", lapply(agg.fun, function(x) {
##                                     y <- cbind(aggregate(object$value,
##                                                          by = list(seqnames = as.factor(seqnames(object)),
##                                                                    key = as.factor(object$key),
##                                                                    start = start(object),
##                                                                    end = end(object),
##                                                                    width = width(object)), FUN = x, ...), x);
##                                     colnames(y) <- c("seqnames", "key", "start", "end", "width", "value", "population");
##                                     y
##                                 }))
##     gr <- GRanges(seqnames = agg.res$seqnames,
##                   ranges = IRanges(agg.res$start, end = agg.res$end), key = agg.res$key,
##                   value = agg.res$value, population = agg.res$population)
##     c(object, gr)
## })
