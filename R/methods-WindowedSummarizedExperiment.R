##' @rdname methods-WindowedSummarizedExperiment
##' @title methods-WindowedSummarizedExperiment
##'
##' @aliases WindowedSummarizedExperiment, methods-WindowedSummarizedExperiment
##'
##' @description Methods that act on objects of class WindowedSummarizedExperiment
##'


##' @rdname window.size
##' @description Get window.size from a WindowedSummarizedExperiment object
##'
##' @param obj WindowedSummarizedExperiment object
##'
setMethod("window.size", "WindowedSummarizedExperiment",
          function(obj) {
    return (rowRanges(obj)@window.size)
})


##' summary.WindowedSummarizedExperiment
##'
##' Summarize the statistics of a WindowedSummarizedExperiment object
##'
##'
##' @param obj WindowedSummarizedExperiment object
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
summary.WindowedSummarizedExperiment <- function(obj, which = NULL, fun = "mean",
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


##' Aggregate WindowedSummarizedExperiment statistics
##'
##' @param x WindowedSummarizedExperiment object
##' @param formula formula
##' @param FUN function to apply
##' @param per.site normalize window scores by window length (kb)
##' @param ... additional arguments passed to aggregate
##'
##' @export
##'
setMethod("aggregate", "WindowedSummarizedExperiment", function(x, formula, FUN, which = NULL, per.site=FALSE, ...) {
    df <- as.data.frame(x, ...)
    res <- aggregate(formula, df, FUN, ...)
    res
})


## override cbind; for combining objects with different ranges
## See https://gist.github.com/PeteHaitch/8993b096cfa7ccd08c13 for discussion
setMethod("cbind", "WindowedSummarizedExperiment",
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
    WindowedSummarizedExperiment(rowRanges = rowRanges(sexp), colData = colData(sexp), assays = assays(sexp))
})


## For plotting; select assay and column and convert to data fram together with rowRanges
setMethod("as.data.frame", "WindowedSummarizedExperiment",
          function(x, row.names = NULL, optional = FALSE,
                   long = FALSE, expand = FALSE,
                   ...) {
    assayData <- as.data.frame(assays(x))
    colnames(assayData)[1:2] <- c("group", "assay")
    if (long)
        assayData <- gather(assayData, select = c(-group, -assay))
    if (expand) {
        i <- match(assayData$key, rownames(colData(x)))
        expand.df <- DataFrame(colData(x)[i,])
        if (ncol(expand.df) == 1)
            colnames(expand.df) <- names(colData(x))
        assayData <- cbind(assayData, as.data.frame(expand.df))
    }
    cbind(rowRanges(x), assayData)
})

setAs("WindowedSummarizedExperiment", "GRanges", function(from) {
    assayData <- as.data.frame(assays(from))
    colnames(assayData)[1:2] <- c("group", "assay")
    gr <- rep(rowRanges(from), nrow(assayData) / length(rowRanges(from)))
    mcols(gr) <- assayData
    if (length(names(mcols(from))) > 0)
        mcols(gr)[names(mcols(from))] <- mcols(from)
    gr
})

##' sites
##'
##'
##' @describeIn sites Get sites for a WindowedSummarizedExperiment object
##'
setMethod("sites", "WindowedSummarizedExperiment",
          function(obj) {sites(rowRanges(obj))})

##' normalize
##'
##' normalize window scores by number of good sites per window
##'
##' @param object WindowedSummarizedExperiment object
##' @param ...
##'
##' @return normalized WindowedSummarizedExperiment object
##' @author Per Unneberg
##'
##' @rdname normalize
##'
##' @importFrom BiocGenerics normalize
##'
setMethod("normalize", "WindowedSummarizedExperiment", function(object, ...) {
    dfl <- DataFrameList(lapply(assayNames(object), function(x){DataFrame(as.matrix(assay(object, x)) / sites(object)) }))
    names(dfl) <- assayNames(object)
    assays(object) <- dfl
    object
})

##' @rdname methods-WindowedSummarizedExperiment
##' @aliases plot
##' @description Plot assay scores versus chromosome position
##' @export
##'
##'
##' @param x A WindowedSummarizedExperiment object
##' @param y optional
##' @param ... arguments to pass to ggplot2
##' @param mapping aesthetic mapping
##' @param size point/line size
##' @param type plot type
##' @param long use long data frame representation of x
##' @param expand expand colData
##' @param wrap wrap plot expression
##' @param strip.position where to align strip
##' @param ncol number of columns in wrap plot
##' @param scales wrap scales
##' @param show.legend whether or not to show legend
##' @param assays what assays to include
##'
##' @import viridis
##'
setMethod("plot", signature=c("WindowedSummarizedExperiment", "missing"),
          function(x, y, ..., mapping = aes(x = pos, y = value), size = .4,
                   type = "p", long = TRUE, expand = TRUE, wrap = NULL,
                   strip.position = "left", ncol = 1, scales = "free_x",
                   show.legend = FALSE, assays = NULL, color.fun = scale_color_viridis,
                   n.levels = 2) {
    type <- match.arg(type, c("p", "l"))
    if (!("index" %in% colnames(as.data.frame(rowRanges(x)))))
        rowRanges(x) <- makeCoordinates(rowRanges(x))

    if (!("colour" %in% names(mapping)))
        rowRanges(x) <- addSeqnamesColor(rowRanges(x), n.levels = n.levels)

    data <- as.data.frame(x, long = long, expand = expand)
    if (!is.null(assays))
        data <- data[data$assay %in% assays, ]
    p <- ggplot(data, mapping = mapping, ...)
    if (!("colour" %in% names(mapping)))
        p <- p + aes(colour = colour)
    if (type == "p") {
        p <- p + geom_point(size = size, show.legend = show.legend)
    } else if (type == "l") {
        p <- p + geom_line(size = size, show.legend = show.legend)
    }
    ## Add color if n.levels
    p <- p + color.fun(n.levels, discrete=TRUE)
    if (is.null(wrap) & long)
        wrap <- . ~ assay
    if (!is.null(wrap))
        p <- p + facet_wrap(wrap, scales = scales,
                            ncol = ncol, strip.position = strip.position)
    p
})

##' @rdname methods-WindowedSummarizedExperiment
##' @aliases plot
##' @description Scatter plot assay scores from different WindowedSummarizedExperiment
##'     objects
##' @export
##' @param x WindowedSummarizedExperiment object
##' @param y WindowedSummarizedExperiment object
##' @param ... parameters to pass to ggplot
##' @param mapping aesthetic mapping
##' @param size plot size
##' @param show.legend show legend
##' @param wrap wrap expression
##' @param scales wrap coordinate scales
##' @param ncol number of wrap columns
##' @param strip.position strip position
##' @param assays which assays to plot
##' @param long use long data frame representations
##' @param expand expand rowRanges columns
##'
setMethod("plot", signature=c("WindowedSummarizedExperiment", "WindowedSummarizedExperiment"),
          function(x, y, ..., mapping = aes(x = x, y = y),
                   size = .5, show.legend = FALSE, wrap = NULL,
                   scales = "fixed", ncol = 1, strip.position = "left",
                   assays = NULL, long = TRUE, expand = TRUE) {
    x <- subsetByOverlaps(x, y)
    y <- subsetByOverlaps(y, x)
    ## Assume identical columns
    z <- x
    assayData <- c(assays(x), assays(y))
    assays(z) <- assayData
    data <- as.data.frame(z, long = long, expand = expand)
    ## Map group to x and y
    i.x <- data$group %in% seq(1, length(assays(x)))
    data$group <- "y"
    data$group[i.x] <- "x"
    if (!is.null(assays))
        data <- data[data$assay %in% assays, ]
    data <- data %>% spread(group, value)
    p <- ggplot(data, mapping = mapping, ...)
    p <- p + geom_point(size = size, show.legend = show.legend)
    if (is.null(wrap) & long)
        wrap <- . ~ assay
    if (!is.null(wrap))
        p <- p + facet_wrap(wrap, scales = scales,
                            ncol = ncol, strip.position = strip.position)
    p
})
