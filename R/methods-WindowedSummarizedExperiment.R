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

##' Aggregate WindowedSummarizedExperiment statistics
##'
##' @param x WindowedSummarizedExperiment object
##' @param formula formula
##' @param FUN function to apply
##' @param ... additional arguments passed to aggregate
##'
##' @export
##'
setMethod("aggregate", "WindowedSummarizedExperiment", function(x, formula, FUN, ...) {
    df <- as.data.frame(x, ...)
    res <- aggregate(formula, df, FUN, ...)
    res
})


##' @rdname cbind
##'
##' @title cbind
##'
##' @description Combine objects with different ranges. For details
##'     and inspiration, see
##'     https://gist.github.com/PeteHaitch/8993b096cfa7ccd08c13 for
##'     discussion. Missing rows are filled with NA values.
##'
##' @param ... Parameters to pass to merge function
##' @param deparse.level See ‘?base::cbind’ for a description of this argument.
##'
##' @return A new WindowedSummarizedExperiment object with merged (union) rowRanges
##' @author Per Unneberg
##'
##'
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



##' @rdname as.data.frame
##' @title as.data.frame
##'
##' @description Function to coerce to a data frame, if possible.
##'
##' @param x The object to coerce
##' @param row.names See ‘?base::as.data.frame’ for a description
##' @param optional See ‘?base::as.data.frame’ for a description
##' @param long Stack assays
##' @param expand Include colData information as extra columns. Only works for long format.
##'
##' @return data.frame
##' @author Per Unneberg
setMethod("as.data.frame", "WindowedSummarizedExperiment",
          function(x, row.names = NULL, optional = FALSE,
                   long = FALSE, expand = FALSE) {
    assaylist <- assayNames(x)[!grepl("^\\..+", assayNames(x))]
    assayData <- as.data.frame(assays(x)[assaylist])

    colnames(assayData)[1:2] <- c("group", "assay")
    if (long) {
        assayData <- gather(assayData, select = c(-group, -assay))
        if (expand) {
            i <- match(assayData$key, rownames(colData(x)))
            expand.df <- DataFrame(colData(x)[i,])
            if (ncol(expand.df) == 1)
                colnames(expand.df) <- names(colData(x))
            assayData <- cbind(assayData, as.data.frame(expand.df))
        }
    }
    cbind(rowRanges(x), assayData)
})


##' sites
##'
##'
##' @describeIn sites Get sites for a WindowedSummarizedExperiment object
##'
setMethod("sites", "WindowedSummarizedExperiment",
          function(obj) {
    if (inherits(rowRanges(obj), "SWindows")) {
        sites(rowRanges(obj))
    } else if (".sites" %in% assayNames(obj)) {
        assay(obj, ".sites")
    } else if (".coverage" %in% assayNames(obj)) {
        y <- DataFrame(apply(assay(obj, ".coverage"), 2, function(x) {x * width(obj)}))
        colnames(y) <- colnames(assay(obj, ".coverage"))
        y
    } else {
        warning("no .coverage or .sites assay for object")
    }
})

##' normalize
##'
##' normalize window scores by number of good sites per window.
##' "hidden" assays whose name begin with a dot are excluded from
##' normalization.
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
    assaylist <- assayNames(object)
    dfl <- DataFrameList(
        lapply(assaylist, function(x) {
            if (grepl("^\\..+", x))
                return(assay(object, x))
            if (is.null(dim(sites(object))))
                DataFrame( as.matrix(assay(object, x)) / (sites(object) %*% t.default(rep(1, ncol(assay(object, x))))))
            else
                DataFrame(as.matrix(assay(object, x)) / as.matrix(sites(object)))
        }))
    names(dfl) <- assaylist
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
