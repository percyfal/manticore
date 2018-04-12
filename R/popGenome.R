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



##'
##' Retrieve genome stats for a PopGenome GENOME instance. Note that
##' this is only intended for use with some slots/statistical
##' operations. The PopGenome function does not work on data that has
##' been saved in rda files. They should be saved with save.session,
##' but this does not seem to work.
##'
##' @export
##' @rdname GStats
##'
##' @examples
##' library(PopGenome)
##' fn <- system.file("extdata", "popgenome.rda", package = "nonmodelr")
##' load(fn)
##' gl <- GENOMEList(scaffolds)
##' sessionInfo()
##' gs <- GStats(gl)
##'
##' @importFrom tidyselect one_of
##' @importFrom utils combn
##'
setMethod("GStats", "GENOMEList",
          function(object,
                   statistics=c("detail"),
                   use.population.names=TRUE,
                   use.region.names=TRUE,
                   quiet=TRUE,
                   ...) {
    if (!requireNamespace("PopGenome", quietly = TRUE))
            stop(paste0("Package \"PopGenome\" needed for this function to work. Please install it."),
                 call. = FALSE)

    statistics <- match.arg(statistics, c("summary", "detail", "neutrality",
                                "fixed", "shared", "diversity",
                                "diversity.between", "F_ST",
                                "F_ST.pairwise", "segregating.sites",
                                "linkage", "sweeps", "recomb"))
    res <- data.frame()
    populations <- names(object[[1]]@populations)
    pairs <- combn(populations, 2, function(x){paste(x, collapse = "/")})
    fun <- list(summary = get.sum.data, detail = get.detail, neutrality = get.neutrality,
                fixed = get.fixed, shared = get.shared,
                diversity = get.diversity, diversity.between = get.diversity.between,
                F_ST = get.F_ST, F_ST.pairwise = get.F_ST.pairwise, linkage = get.linkage,
                sweeps = get.sweeps, recomb = get.recomb, segregating.sites = get.segregating.sites)

    f <- fun[[statistics]]
    gather.key = "key"
    for (name in names(object)) {
        ranges <- object[[name]]@region.names
        if (statistics %in% c("detail", "neutrality", "diversity", "linkage", "sweeps", "recomb")) {
            ## Result is matrix with populations as row names
            tmp <- f(object[[name]])
            if (use.population.names) {
                rownames(tmp) <- populations
            }
            tmp <- do.call("rbind", lapply(rownames(tmp), function(x){data.frame(population = x, ranges = rownames(tmp[x,][[1]]), tmp[x,][[1]])}))
            tmp$seqnames <- name
            rownames(tmp) <- NULL
            gather.exclude <- c("population", "ranges", "seqnames")
        } else if (statistics %in% c("fixed", "shared", "diversity.between")) {
            ## Result is data frame with population pairs in columns
            tmp <- as.data.frame(f(object[[name]], ...))
            if (use.population.names) {
                colnames(tmp) <- pairs
            }
            tmp$seqnames <- name
            tmp$ranges <- ranges
            tmp$key <- statistics
            gather.key <- "population"
            gather.exclude <- c("ranges", "seqnames", "key")
        } else if (statistics %in% c("F_ST.pairwise")) {
            ## Result is matrix with statistics as row names
            tmp <- f(object[[name]])
            col.names <- colnames(tmp[[1]])
            tmp <- do.call(
                "rbind", lapply(rownames(tmp),
                                function(x) {
                             data.frame(key = x, ranges = ranges, tmp[x, ])}))
            if (use.population.names) {
                colnames(tmp)[3:dim(tmp)[2]] <- pairs
            } else {
                colnames(tmp)[3:dim(tmp)[2]] <- col.names
            }
            tmp$seqnames <- name
            gather.key <- "population"
            gather.exclude <- c("key", "ranges", "seqnames")
        } else if (statistics %in% c("summary")) {
            ## Result is matrix with statistics in columns, region in rows
            if (!quiet) message("Analyzing ", statistics, " data")
            tmp <- as.data.frame(f(object[[name]]))
            tmp$seqnames <- name
            tmp$ranges <- ranges
            tmp$population <- "__all__"
            gather.exclude <- c("seqnames", "ranges", "population")
        } else if (statistics %in% c("F_ST")) {
            ## Result is matrix with statistics in column, regions in rows
            tmp <- as.data.frame(f(object[[name]]))
            tmp$ranges <- ranges
            tmp$seqnames <- name
            tmp$population <- "__all__"
            gather.exclude <- c("ranges", "seqnames", "population")
        } else if (statistics %in% c("segregating.sites")) {
            ## Result is matrix with populations in columns
            tmp <- as.data.frame(f(object[[name]]))
            if (use.population.names) {
                colnames(tmp) <- populations
            }
            tmp$ranges <- ranges
            tmp$seqnames <- name
            tmp$key <- statistics
            gather.key <- "population"
            gather.exclude <- c("ranges", "seqnames", "key")
        } else {
            stop("shouldn't end up here")
        }
        if (!quiet) message("Adding data for ", name)
        res <- rbind(res, tmp)
    }

    ## Fix factors
    res$seqnames <- factor(res$seqnames, levels = unique(res$seqnames))
    res$ranges <- factor(res$ranges, levels = unique(res$ranges))
    res <- res %>% tidyr::gather(gather.key, "value", -one_of(gather.exclude))
    colnames(res)[which(colnames(res) == "gather.key")] <- gather.key
    res$population <- factor(res$population, levels = unique(res$population))
    res$key <- factor(res$key, levels = unique(res$key))

    data <- res %>% tidyr::unite("population_statistics", population, key) %>%
        tidyr::spread("population_statistics", value)
    ## Parse the start/stop positions; we assume the format "name start - end :"
    data.ranges <- strsplit(as.character(data$ranges), " ")
    if (all(lapply(data.ranges, length) == 5)) {
        start <- unlist(lapply(data.ranges, function(x) {as.integer(x[[2]])}))
        end <- unlist(lapply(data.ranges, function(x) {as.integer(x[[4]])}))
    } else {
        start <- as.integer(rep(1, length(data.ranges)))
        end <- as.integer(unlist(lapply(data$seqnames, function(x) {object[[as.character(x)]]@n.sites})))
    }
    .values <- subset(data, select = -c(ranges, seqnames))
    .ranges <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(data$seqnames, rep(1, length(data$seqnames))),
                                      ranges = IRanges::IRanges(start, end = end),
                                      feature_id = paste(data$seqnames, start, end, ":", sep = " "))
    .ranges$feature_id <- factor(.ranges$feature_id, levels=unique(.ranges$feature_id))
    rse <- SummarizedExperiment(assays = list(data = as.matrix(.values)),
                                rowRanges = .ranges)
    .cdata <- S4Vectors::DataFrame(
        expand.grid(
            statistic = sort(levels(res$key), decreasing = FALSE),
            population = sort(levels(res$population), decreasing = FALSE)),
        row.names = colnames(rse))
    colData(rse) <- .cdata
    gs <- new("GStats", rse, statistics = statistics, application = "PopGenome")
    return (gs)
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


##' genomewide.stats
##'
##' @export
##' @rdname genomewide.stats
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
