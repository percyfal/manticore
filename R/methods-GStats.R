##' summary.GStats
##'
##' Summarize the statistics of a GStats object by regions or
##' populations
##'
##' @param gs GStats object
##' @param per.site calculate per site values (kb)
##' @param fun summary function
##' @param per.region normalize score by window size
##' @param ... additional arguments
##' @return a data frame
##' @author Per Unneberg
##' @export
##'
summary.GStats <- function(gs, per.site=FALSE, fun="mean", per.region=FALSE, ...) {
    fun <- match.arg(fun, c("mean", "sd", "median", "var"))
    df <- assay(gs)
    if (per.site) {
        sites <- rowRanges(gs)$sites
        sites[sites == 0] <- NA
        df <- df / sites
    }
    df[is.infinite(df)] <- NA
    df[is.nan(as.matrix(df))] <- NA
    if (!per.region) {
        m <- matrix(apply(df, 2, fun, ...), nrow = length(levels(gs$population)), byrow = TRUE)
        ret <- as.data.frame(m)
        rownames(ret) <- levels(gs$population)
        colnames(ret) <- levels(gs$statistic)
    } else {
        m <- do.call("rbind", by(t(df), gs$statistic, function(x) {apply(x, 2, fun, ...)}))
        ret <- as.data.frame(t(m), row.names = as.character(SummarizedExperiment::rowData(gs)$feature_id))
    }
    ret
}


##' Aggregate GStats statistics
##'
##' @param x GStats object
##' @param formula formula
##' @param FUN function to apply
##' @param per.site normalize window scores by window length (kb)
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
##' @param per.site calculate per site statistics (divide by width) per kb
##' @return GRanges instance
##' @author Per Unneberg
##'
##' @rdname asGRanges
##' @export
##'
##' @importFrom GenomicRanges GRanges
##' @importFrom IRanges IRanges
##' @importFrom S4Vectors values
##'
setMethod("asGRanges", "GStats", function(x, long=TRUE, per.site=FALSE) {
    gr <- rowRanges(x)
    y <- assay(x)
    if (per.site) {
        sites <- rowRanges(x)$sites
        y <- y / sites
    }
    y[is.infinite(y)] <- NA
    y[is.nan(y)] <- NA
    y <- as.data.frame(y)
    colnames(y) <- make.names(rownames(colData(x)))
    values(gr) <- cbind(as.data.frame(values(gr)), y)
    if (long) {
        exclude <- names(elementMetadata(rowRanges(x)))
        df <- tidyr::gather(as.data.frame(values(gr)), statistic, value, - exclude)
        i <- match(df$statistic, make.names(rownames(SummarizedExperiment::colData(x))))
        if (any(is.na(i)))
            stop("names in long data frame do not match those of rownames(colData(x)); check for special characters in population names")
        n.rep <- nrow(df) / length(gr)
        gr <- GRanges(
            ranges = IRanges(start = rep(start(gr), n.rep),
                             end = rep(end(gr), n.rep)),
            seqnames = rep(seqnames(gr), n.rep),
            strand = rep(strand(gr), n.rep),
            seqinfo = seqinfo(gr))
        values(gr) <- df
        gr$statistic <- x$statistic[i]
        gr$population <- x$population[i]
    }
    gr
})


##' Create GStats from GRanges object
##'
##' @export
##'
##' @importFrom S4Vectors Rle
setMethod("GStats", signature(object = "GRanges"),
          function(object,
                   population,
                   statistics,
                   application=NA,
                   ...) {
    .values <- values(object)
    .ranges <- GRanges(seqnames=seqnames(object), ranges=ranges(object), feature_id = paste(seqnames(object), start(object), end(object), ":", sep=" "), sites=0)
    .ranges$feature_id <- factor(.ranges$feature_id, levels = unique(.ranges$feature_id))
    rse <- SummarizedExperiment(assays = list(data = as.matrix(.values)),
                                rowRanges = .ranges)
    .cdata <- S4Vectors::DataFrame(
                             expand.grid(
                                 statistic = sort(levels(factor(statistics)), decreasing = FALSE),
                                 population = sort(levels(factor(population)), decreasing = FALSE)),
                             row.names = colnames(rse))
    colData(rse) <- .cdata
    gs <- new("GStats", rse, statistics = statistics, application = application)
    return (gs)
})
