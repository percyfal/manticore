##' @rdname popGenome
##' @title popGenome
##'
##' @description Methods and utilities for dealing with popGenome
##'     results
##'


##' GENOMEList
##'
##' Convert list to GENOMEList
##'
##' @export
##' @rdname GENOMEList
##'
##' @importFrom methods new
##'
setMethod("GENOMEList", "list",
          function(obj) new("GENOMEList", listData = obj))

##' GRanges_OR_NULL
##'
##' Create class union of GRanges or NULL
##'
##' @export
##' @rdname GRanges_OR_missing
##'
##' @importFrom methods setClassUnion
##' @importFrom GenomicRanges GRanges
##'
setClassUnion("GRanges_OR_missing", c("GRanges", "missing"))


##' @rdname Seqinfo_OR_missing
##' @title Seqinfo_OR_missing
##'
##' @description Class union of Seqinfo or missing
##'
##' @export
##'
##' @importFrom methods setClassUnion
##' @import GenomeInfoDb
##'
setClassUnion("Seqinfo_OR_missing", c("Seqinfo", "missing"))



## Interntal functions
.population.data.frame <- function(x) {
    y <- do.call(rbind, x)
    n <- nrow(y) / length(rownames(x))
    cbind(population = rep(rownames(x), each = n), seqnames = rownames(y), y)
}

.make.population.assay.list <- function(x) {
    y <- do.call(rbind, lapply(x, .population.data.frame))
    z <- as.data.frame(y) %>% gather("key", "value", -one_of(c("seqnames", "population"))) %>% tidyr::spread("population", value)
    lapply(split(z, z$key), function(x) {y <- as.data.frame(x)[, !(colnames(x) %in% exclude)]; rownames(y) <- NULL; y})
}


.pairs.data.frame <- function(x) {
    y <- do.call(rbind, x)
    n <- nrow(y) / length(rownames(x))
    cbind(key = rep(rownames(x), each = n), y)
}

.make.pairs.assay.list <- function(x, stat) {
    if (stat == "F_ST.pairwise")
        z <- as.data.frame(do.call("rbind", lapply(x, .pairs.data.frame)))
    else
        z <- as.data.frame(cbind(key = stat, do.call(rbind, x)))
    lapply(split(z, z$key), function(x) {y <- as.data.frame(x)[, !(colnames(x) %in% exclude)]; rownames(y) <- NULL; y})
}


##' @rdname PopGenome
##'@title Popgenome
##' @description Function for reading PopGenome results
##'
##' @param object an R object
##'
setGeneric("PopGenome",
           function(object, seqinfo, ...) standardGeneric("PopGenome"))


##' @rdname PopGenome
##'
##'
##' @description Process a GENOMEList object and return a list
##'     containing WindowedSummarizedExperiment objects
##'
##' @param object GENOMEList
##' @param statistics Which statistics to parse
##' @param quiet don't print info messages
##' @param window.size Set window size if unable to infer
##'
##' @examples
##' library(PopGenome)
##' fn <- system.file("extdata", "popgenome.rda", package = "manticore")
##' load(fn)
##' wselist <- PopGenome(gl)
##'
##' @return list
##' @author Per Unneberg
##'
setMethod("PopGenome", signature("GENOMEList", "Seqinfo_OR_missing"),
          function(object,
                   seqinfo = NULL,
                   statistics = c("detail", "neutrality", "fixed", "shared", "diversity",
                                "diversity.between", "F_ST", "F_ST.pairwise"),
                   quiet = TRUE,
                   window.size = NULL) {
    if (!requireNamespace("PopGenome", quietly = TRUE))
            stop(paste0("Package \"PopGenome\" needed for this function to work. Please install it."),
                 call. = FALSE)
    statistics <- match.arg(statistics, c("detail", "neutrality",
                                          "fixed", "shared",
                                          "diversity",
                                          "diversity.between", "F_ST",
                                          "F_ST.pairwise",
                                          "linkage", "sweeps",
                                          "recomb"),
                            several.ok = TRUE)
    population.statistics <- c("detail", "neutrality", "diversity", "linkage", "sweeps", "recomb")
    pairs.statistics <- c("fixed", "shared", "diversity.between", "F_ST.pairwise")
    fun <- list(detail = get.detail, neutrality = get.neutrality,
                fixed = get.fixed, shared = get.shared,
                diversity = get.diversity, diversity.between = get.diversity.between,
                F_ST = get.F_ST, F_ST.pairwise = get.F_ST.pairwise, linkage = get.linkage,
                sweeps = get.sweeps, recomb = get.recomb)

    assayData <- S4Vectors::SimpleList()
    colData <- S4Vectors::DataFrame(populations = names(object[[1]]@populations))

    ## Generate WSE for populations data
    rowRanges.df <- DataFrame()
    for (name in names(object)) {
        obj <- object[[name]]

        ## Find out if slide data or not
        if (length(obj@SLIDE.POS) == 0) {
            if (!quiet) message("collecting summary data")
            rowRanges.df <- rbind(rowRanges.df, DataFrame(seqnames = name, start = 1, end = get.n.sites(obj)))
        } else {
            ranges <- as.data.frame(do.call("rbind", strsplit(obj@region.names, " ")),
                                    stringsAsFactors = FALSE)[c(1,2,4)]
            colnames(ranges) <- c("seqnames", "start", "end")
            rowRanges.df <- rbind(rowRanges.df, DataFrame(ranges))
        }
    }

    if (!quiet) message("collecting segregating sites")
    segregating.sites <- as.data.frame(do.call(rbind, lapply(object, get.segregating.sites)))
    assayData$segregating.sites <- segregating.sites

    for (stat in statistics) {
        if (!(stat %in% population.statistics))
            next
        if (!quiet) message("collecting ", stat, " population results")
       assayData <- c(assayData, .make.population.assay.list(lapply(object, fun[[stat]])))
    }
    if (is.null(window.size)) {
        window.size <- as.integer(ranges$end[1]) - as.integer(ranges$start[1]) + 1
        message("window.size parameter undefined; inferring window size to ", window.size, " from data")
    }
    rowRanges <- Windows(seqnames = rowRanges.df$seqnames, ranges = IRanges(start = as.integer(rowRanges.df$start),
                                                                            end = as.integer(rowRanges.df$end)),
                         window.size = window.size)

    if (!is.null(seqinfo))
        seqinfo(rowRanges) <- seqinfo
    wse.population <- WindowedSummarizedExperiment(assays = assayData, rowRanges = rowRanges, colData = colData)

    ## Get population pairs
    poppairs <- combn(names(object[[1]]@populations), 2, function(x){paste(x, collapse = "/")})
    colData <- S4Vectors::DataFrame(population.pairs = poppairs)
    rownames(colData) <- combn(seq(length(poppairs)), 2, function(x) {paste0("pop", x[1], "/", "pop", x[2])})
    assayData <- S4Vectors::SimpleList()
    for (stat in statistics) {
        if (!(stat %in% pairs.statistics))
            next
        if (!quiet) message("collecting ", stat, " pairwise results")
        assayData <- c(assayData, .make.pairs.assay.list(lapply(object, fun[[stat]]), stat))
    }
    wse.pair <- WindowedSummarizedExperiment(assays = assayData, rowRanges = rowRanges, colData = colData)


    colData <- S4Vectors::DataFrame(all = "all")
    assayData <- S4Vectors::SimpleList()
    if ("F_ST" %in% statistics) {
        message("Class of object: ", class(object))
        data <- as.list(as.data.frame(do.call(rbind, lapply(object, fun[["F_ST"]]))))
        data.list <- lapply(names(data), function(x){y <- DataFrame(all = data[[x]])})
        names(data.list) <- names(data)
        assayData <- c(assayData, data.list)
    }
    wse.all <- WindowedSummarizedExperiment(assays = assayData, rowRanges = rowRanges, colData = colData)
    list(population = wse.population,
         pair = wse.pair,
         all = wse.all)
})


## No need to make these generic; they are not exported
get.fixed <- function(object, ...) {
    return (slot(object, "n.fixed.sites"))
}

get.shared <- function(object, ...) {
    return (slot(object, "n.shared.sites"))
}

get.diversity.between <- function(object, which="nuc", ...) {
    which <- match.arg(which, c("nuc", "hap"))
    slots <- list(nuc = "nuc.diversity.between", hap = "hap.diversity.between")
    df <- slot(object, slots[[which]])
    if (is.null(colnames(df))) {
        df <- t(df)
    }
    ncol <- dim(df)[2]
    col.names <- colnames(df)
    df <- data.frame(matrix(unlist(df), ncol = ncol,
                            byrow = TRUE), stringsAsFactors = FALSE)
    colnames(df) <- col.names
    return (df)
}

get.F_ST.pairwise <- function(object, which, ...) {
    df <- get.F_ST(object, pairwise = TRUE)
    return (df)
}

get.segregating.sites <- function(object, ...) {
    return (object@n.segregating.sites)
}

get.n.sites <- function(object) {
    return (get.sum.data(object)[1])
}



##' @rdname genomewide.stats
##'
##' @export
##'
setGeneric("genomewide.stats", function(object, which=character(0), biallelic.structure=TRUE,
                                        pi=TRUE, ...) standardGeneric("genomewide.stats"))
##' genomewide.stats
##'
##' Calculate genome wide statistics
##'
##' @param object GENOME object on which to perform calculations
##' @param which statistics to calculate
##' @param biallelic.structure calculate biallelic structure
##' @param pi Nei's calculation of pi
##' @param ... extra arguments
##' @return Updated object
##' @author Per Unneberg
##' @export
##' @rdname genomewide.stats
##'
setMethod("genomewide.stats", "GENOME", function(object,
                                                 which=c("detail", "neutrality", "F_ST",
                                                         "diversity", "diversity.between",
                                                         "fixed.shared"),
                                                 biallelic.structure=TRUE, pi=TRUE,
                                                 ...) {
    which <- match.arg(which, c("detail", "neutrality", "F_ST", "diversity", "diversity.between", "fixed.shared", "linkage", "R2", "recomb", "sweeps"), several.ok = TRUE)

    functions <- list(detail = detail.stats, neutrality = neutrality.stats, F_ST = F_ST.stats,
                      diversity = diversity.stats, diversity.between = diversity.stats.between,
                      fixed.shared = calc.fixed.shared, linkage = linkage.stats, R2 = calc.R2,
                      recomb = recomb.stats, sweeps = sweeps.stats)

    call.genomewide.stats.function <- function(object, fn, ...) {
        message(paste0("\napplying analysis ", fn))
        object <- tryCatch({
            object <- functions[[fn]](object, ...)
        }, error = function(err) {
            message(paste0("\nFailed to apply analysis ", fn))
            message(err)
            return(object)
        })
        object
    }
    if ("detail" %in% which) {
        object <- call.genomewide.stats.function(object, "detail", biallelic.structure = biallelic.structure, ...)
    }
    if ("neutrality" %in% which) {
        object <- call.genomewide.stats.function(object, "neutrality", ...)
    }
    if ("F_ST" %in% which) {
        object <- call.genomewide.stats.function(object, "F_ST", ...)
    }
    if ("diversity" %in% which) {
        object <- call.genomewide.stats.function(object, "diversity", pi = pi, ...)
    }
    if ("diversity.between" %in% which) {
        object <- call.genomewide.stats.function(object, "diversity.between", ...)
    }
    if ("fixed.shared" %in% which) {
        object <- call.genomewide.stats.function(object, "fixed.shared", ...)
    }
    if ("linkage" %in% which) {
        object <- call.genomewide.stats.function(object, "linkage", ...)
    }
    if ("R2" %in% which) {
        object <- call.genomewide.stats.function(object, "R2", ...)
    }
    if ("recomb" %in% which) {
        object <- call.genomewide.stats.function(object, "recomb", ...)
    }
    if ("sweeps" %in% which) {
        object <- call.genomewide.stats.function(object, "sweeps", ...)
    }
    return (object)
})

##' plot_region_stats
##'
##' Plot region stats
##' @title plot_region_stats
##' @param object object
##' @param x.var a name identifying a variable to plot values against
##' @param wrap wrap plot
##' @param agg.fun aggregation functions
##' @param text.size text size
##' @param which which plots to include. Choose from population levels and aggregation function
##' @param label.outlier label outlier points
##' @param label.format use short (default) or long format for window labels
##' @param n.se number of standard errors to use as cutoff for identifying outliers, based on confidence interval level
##' @param label.hjust horizontal adjustment for labels
##' @param label.vjust vertical adjustment for labels
##' @param level level of confidence interval to use
##' @param method method to use for fitting
##' @param scales scales for facet plots
##' @param ... additional arguments to pass to geom_point()
##' @return ggplot
##' @author Per Unneberg
##' @export
setGeneric("plot_region_stats", function(object, x.var="width", wrap=TRUE, agg.fun=c("sum", "mean"), text.size=12, which=NULL, label.outlier=FALSE, label.format="short", n.se=5, label.hjust=1, label.vjust=1, level=0.99, method="loess", scales="free_y", ...) standardGeneric("plot_region_stats"))
##' @describeIn plot_region_stats Plot statistics versus region for GRanges instance
setMethod("plot_region_stats", "GRanges", function(object, x.var, wrap, agg.fun, text.size, which, label.outlier, n.se, label.hjust, label.vjust, level, method, scales, ...) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop(paste(
            "Package \"GenomicRanges\" needed for this function to work.",
            "Please install it.", sep = ""),
            call. = FALSE)
    }
    x.var <- match.arg(x.var, c("width", "seqlengths", setdiff(colnames(values(object)), c("key", "value", "population"))))
    label.format <- match.arg(label.format, c("long", "short"))
    if (x.var == "seqlengths") {
        if (any(is.na(seqlengths(object)))) stop("no seqlengths defined for object; must be set for this function to work")
        values(object)[[x.var]] <- seqlengths(object)
    }
    which.options <- levels(as.factor(object$population))
    if (is.null(which)) {
        which <- which.options
    } else {
        which <- match.arg(which, which.options)
    }
    if (label.outlier & !("outlier" %in% colnames(values(object)))) {
        message("no outlier column in values data frame; calling function identify_outliers")
        conf <- identify_outliers(as.data.frame(object), paste("value ~ ", x.var), key="population", method=method, level=level)
        if (label.format=="long") {
            label <- paste(seqnames(object), start(object), end(object), sep="-")
        } else {
            label <- seqnames(object)
        }
        object$outlier <- ifelse(abs(conf$residuals/conf$se.fit) > n.se, label, "")
    }
    p <- ggplot(subset(as.data.frame(object), population %in% which), aes_string(x = x.var, y = "value")) + geom_point(...) + theme(text = element_text(size = text.size))
    if (label.outlier) p <- p + geom_text(aes(label = outlier), hjust = label.hjust, vjust = label.vjust)
    if (length(which) == 1) wrap <- FALSE
    if (wrap) p <- p + facet_wrap(~ population, scales = scales)
    p
})
