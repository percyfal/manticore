##' @describeIn GENOMEList Convert list to GENOMEList
##' @export
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
##' # Read a vcf and generate two genome objects
##' vcf_file <- system.file("extdata", "medium.call.biallelic.vcf.gz", package="nonmodelr")
##' scaffold1 <- readVCF(vcf_file, 1000, frompos=1, topos=1000000, tid="scaffold1")
##' scaffold13 <- readVCF(vcf_file, 1000, frompos=1, topos=1000000, tid="scaffold13")
##' gl <- GENOMEList(list(scaffold1=scaffold1, scaffold13=scaffold13))
##' gs <- GStats(gl)
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
    res <- res %>% gather(gather.key, "value", -one_of(gather.exclude))
    colnames(res)[which(colnames(res) == "gather.key")] <- gather.key
    res$population <- factor(res$population, levels = unique(res$population))
    res$key <- factor(res$key, levels = unique(res$key))

    data <- res %>% unite("population_statistics", population, key) %>%
        spread("population_statistics", value)
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
    .ranges <- GRanges(seqnames = Rle(data$seqnames, rep(1, length(data$seqnames))),
                       ranges = IRanges(start, end = end),
                       feature_id = paste(data$seqnames, start, end, ":", sep=" "))
    rse <- SummarizedExperiment(assays = list(data = as.matrix(.values)),
                                rowRanges = .ranges)
    .cdata <- DataFrame(
        expand.grid(population = levels(res$population),
                    statistic = levels(res$key)),
        row.names = colnames(rse))
    colData(rse) <- .cdata
    gs <- new("GStats", rse, statistics = statistics)
    return(gs)
})



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
##' Calculate genome wide statistics
##' @title genomewide.stats
##' @param object object on which to perform calculations
##' @param which statistics to calculate
##' @param biallelic.structure calculate biallelic structure
##' @param pi Nei's calculation of pi
##' @param ... extra arguments
##' @return Updated object
##' @author Per Unneberg
##' @export
setGeneric("genomewide.stats", function(object, which=c("detail", "neutrality", "F_ST", "diversity", "diversity.between", "fixed.shared", "linkage"), biallelic.structure=TRUE, pi=TRUE, ...) standardGeneric("genomewide.stats"))
##' @describeIn genomewide.stats Calculate genome wide statistics for a GENOME object
setMethod("genomewide.stats", "GENOME", function(object, which, ...) {
    which <- match.arg(which, c("detail", "neutrality", "F_ST", "diversity", "diversity.between", "fixed.shared", "linkage", "R2", "recomb", "sweeps"), several.ok = TRUE)
    if ("detail" %in% which) {
        object <- detail.stats(object, biallelic.structure = biallelic.structure, ...)
    }
    if ("neutrality" %in% which) {
        object <- neutrality.stats(object, ...)
    }
    if ("F_ST" %in% which) {
        object <- F_ST.stats(object, ...)
    }
    if ("diversity" %in% which) {
        object <- diversity.stats(object, pi = pi, ...)
    }
    if ("diversity.between" %in% which) {
        object <- diversity.stats.between(object, ...)
    }
    if ("fixed.shared" %in% which) {
        object <- calc.fixed.shared(object, ...)
    }
    if ("linkage" %in% which) {
        object <- linkage.stats(object, ...)
    }
    if ("R2" %in% which) {
        object <- calc.R2(object, ...)
    }
    if ("recomb" %in% which) {
        object <- recomb.stats(object, ...)
    }
    if ("sweeps" %in% which) {
        object <- sweeps.stats(object, ...)
    }
    return (object)
})
##' summary.GRanges
##'
##' Summarize the statistics of a GRanges object
##' @title summary.GRanges
##' @param gr GRanges object
##' @param per.site calculate per site values
##' @param ... additional arguments
##' @return a data frame
##' @author Per Unneberg
##' @export
summary.GRanges <- function(gr, per.site=TRUE, fun="mean", ...) {
    fun <- match.arg(fun, c("mean", "sd"))
    df <- as.data.frame(gr)
    if (per.site) df$value <- df$value / width(gr)
    tab <- tapply(df$value, list(df$population, df$key), fun)
    as.data.frame(tab)
}


##' plot.pg
##'
##' Generic plotting function for PopGenome results
##' @title plot.pg
##' @param data long format from GStats function
##' @param x variable to map to x aestethic
##' @param y variable to map to y aestethic
##' @param colour colours to use
##' @param colour.var variable to map to colour aestethic
##' @param wrap wrap plots
##' @param wrap.formula wrap formula
##' @param wrap.ncol number of columns in facet wrap
##' @param plot.type plot type, either point or line
##' @param x.lab x label
##' @param y.lab y label
##' @param main plot title
##' @param compact.facet compact facet representation
##' @param strip.position strip position
##' @param scales scales
##' @param size plot size
##' @param hide.legend whether or not to hide legend
##' @param hide.xaxis hide x axis tick marks and labels
##' @param grid include grid lines
##' @param text.size text size
##' @param text.x.angle x text angle
##' @param text.x.hjust x text horizontal justification
##' @param ... extra arguments
##' @param which statistic to plot
##' @return ggplot
##' @author Per Unneberg
##' @export
plot.pg <- function(data, x="ranges", y="value",
                    colour=brewer.pal(3, "Dark2"), colour.var="seqnames",
                    wrap=TRUE, wrap.formula="key ~ population",
                    wrap.ncol=1, plot.type="point", x.lab="window",
                    y.lab=NULL, main=NULL,
                    compact.facet=TRUE, strip.position="right",
                    scales="free_y", size=1, hide.legend=TRUE,
                    hide.xaxis=TRUE, grid=FALSE,
                    text.size=14, text.x.angle=45, text.x.hjust=1,
                    ...) {
    plot.type <- match.arg(plot.type, c("point", "line"))
    data[[colour.var]] <- factor(data[[colour.var]], levels=unique(data[[colour.var]]))
    ## Make sure factor is ordered according to order of occurrence if not numeric
    if (!is.numeric(data[[x]])) {
        if (!is.factor(data[[x]])) data[[x]] <- factor(data[[x]], levels=unique(data[[x]]))
    }

    p <- ggplot(data, aes_string(x = x, y = y, colour = colour.var))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (plot.type == "point") {
        p <- p + geom_point(size = size)
    } else if (plot.type == "line") {
        p <- p + geom_line(size = size)
    }
    if (!is.null(x.lab)) p <- p + xlab(x.lab)
    if (!is.null(y.lab)) p <- p + ylab(y.lab)
    if (!is.null(main)) p <- p + ggtitle(main)
    if (compact.facet) {
        p <- p + theme(panel.spacing = unit(0, "lines"))
    }
    if (hide.legend) p <- p + theme(legend.position = "none")
    if (!grid) p <- p + theme(panel.grid = element_blank())
    if (hide.xaxis) {
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
        p <- p + theme(axis.text.x = element_text(angle = text.x.angle, hjust = text.x.hjust))
    }
    p <- p + theme(text = element_text(size = text.size))
    # Colouring setup
    nlevels <- length(levels(data[[colour.var]]))
    nc <- length(colour)
    n <- floor(nlevels / nc)
    nmod <- nlevels %% nc
    if (nmod == 0) {
        p <- p + scale_colour_manual(values = rep(colour, n))
    } else {
        p <- p + scale_colour_manual(values = c(rep(colour, n), colour[1:nmod]))
    }
    p
}
##' @title plot.GRanges
##' @describeIn plot.pg Make plot of GRanges object
##' @param per.site normalize statistics to per-site values
##' @export
plot.GRanges <- function(data, which=levels(factor(data$key)), per.site=TRUE, y="value", ...) {
    which <- match.arg(which, levels(factor(data$key)), several.ok = TRUE)
    df <- cbind(seqnames = as.character(seqnames(data)), as.data.frame(values(data)))
    df$ranges <- paste(as.character(seqnames(data)), start(data), "-", end(data), ":", sep = " ")
    if (per.site) df[[y]] <- df[[y]] / width(data)
    plot.pg(subset(df, key %in% which), y = y, ...)
}
##' @title plot.pg.summary
##' @describeIn plot.pg Make plot of summary
##' @export
plot.pg.summary <- function(data, x="seqnames", y="value", wrap.formula="~ key", which=c("n.sites", "trans.transv.ratio"), x.lab="scaffold", hide.xaxis=FALSE, ...) {
    which <- match.arg(which, c("n.sites", "n.biallelic.sites", "n.gaps", "n.unknowns", "n.valid.sites", "n.polyallelic.sites", "trans.transv.ratio"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), x = x, y = y, wrap.formula = wrap.formula, x.lab = x.lab, hide.xaxis = hide.xaxis, ...)
}
##' @title plot.pg.detail
##' @describeIn plot.pg Make plot of details
##' @export
plot.pg.detail <- function(data, which=c("MDG1", "MDG2", "MDSD"), ...) {
    which = match.arg(which, c("MDG1", "MDG2", "MDSD"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.neutrality
##' @describeIn plot.pg Make plot of neutrality statistics
##' @export
plot.pg.neutrality <- function(data, which=c("Tajima.D", "Fu.Li.F", "Fu.Li.F"), ...) {
    which = match.arg(which, c("Tajima.D", "n.segregating.sites", "Rozas.R_2", "Fu.Li.F", "Fu.Li.D", "Fu.F_S", "Fay.Wu.H", "Zeng.E", "Strobeck.S"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.fixed
##' @describeIn plot.pg Make plot of fixed sites
##' @export
plot.pg.fixed <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.shared
##' @describeIn plot.pg Make plot of shared sites
##' @export
plot.pg.shared <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.diversity
##' @describeIn plot.pg Make plot of diversity
##' @export
plot.pg.diversity <- function(data, which=c("nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"),  ...) {
    which <- match.arg(which, c("hap.diversity.within", "hap.F_ST.vs.all", "nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.diversity.between
##' @describeIn plot.pg Make plot of between population diversity
##' @export
plot.pg.diversity.between <- function(data, wrap.formula="~ population", colour="name", ...) {
    plot.pg(data, wrap.formula = wrap.formula, colour = colour, ...)
}
##' @title plot.pg.F_ST
##' @describeIn plot.pg Make plot of F_ST
##' @export
plot.pg.F_ST <- function(data, wrap.formula="~ key", which=c("nucleotide.F_ST", "Nei.G_ST"), ...) {
    which <- match.arg(which, c("haplotype.F_ST", "nucleotide.F_ST", "Nei.G_ST", "Hudson.G_ST", "Hudson.H_ST", "Hudson.K_ST"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.F_ST.pairwise
##' @describeIn plot.pg Make plot of pairwise F_ST
##' @export
plot.pg.F_ST.pairwise <- function(data, ...) {
    plot.pg(data, ...)
}
##' @title plot.pg.segregating.sites
##' @describeIn plot.pg Make plot of segregating sites
##' @export
plot.pg.segregating.sites <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}
##' aggregate_region_stats
##'
##' Aggregate statistics over regions, where the regions typically are
##' scaffolds or windows. Aggregation is performed over the list of
##' components seqnames, key (i.e. statistic), region start, region
##' end, and region width
##' @title aggregate_region_stats
##' @param object object
##' @param agg.fun aggregation functions
##' @param ... additional arguments
##' @return GRanges with augmented values
##' @author Per Unneberg
##' @export
setGeneric("aggregate_region_stats", function(object, agg.fun=c("sum", "mean"), ...) standardGeneric("aggregate_region_stats"))
##' destribeIn aggregate_region_stats Aggregate statistics for GRanges instance
setMethod("aggregate_region_stats", "GRanges", function(object, agg.fun=c("sum", "mean"), ...) {
    agg.res <- do.call("rbind", lapply(agg.fun, function(x) {
                                    y <- cbind(aggregate(object$value,
                                                         by = list(seqnames = as.factor(seqnames(object)),
                                                                   key = as.factor(object$key),
                                                                   start = start(object),
                                                                   end = end(object),
                                                                   width = width(object)), FUN = x, ...), x);
                                    colnames(y) <- c("seqnames", "key", "start", "end", "width", "value", "population");
                                    y
                                }))
    gr <- GRanges(seqnames = agg.res$seqnames,
                  ranges = IRanges(agg.res$start, end = agg.res$end), key = agg.res$key,
                  value = agg.res$value, population = agg.res$population)
    c(object, gr)
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

##' boxplot.pg
##'
##' Make a box/violin plot of data
##' @title boxplot.pg
##' @param formula a formula, such as y ~ grp
##' @param data long format from GStats function
##' @param colour colours to use
##' @param colour.var variable to map to colour aestethic
##' @param wrap wrap plots
##' @param wrap.formula wrap formula
##' @param wrap.ncol number of columns in facet wrap
##' @param plot.type plot type, either point or line
##' @param x.lab x label
##' @param y.lab y label
##' @param main plot title
##' @param compact.facet compact facet representation
##' @param strip.position strip position
##' @param scales scales
##' @param hide.legend whether or not to hide legend
##' @param grid include grid lines
##' @param text.size text size
##' @param text.x.angle angle of x tick labels
##' @param text.x.hjust horizontal adjustment of x tick labels
##' @param ... arguments passed on to facet_wrap
##' @param size plot size
##' @param which statistic to plot
##' @importFrom graphics boxplot
##' @return ggplot object
##' @author Per Unneberg
##' @export
boxplot.pg <- function(formula = "value ~ population", data=NULL,
                       colour=brewer.pal(3, "Dark2"), colour.var=NULL,
                       wrap=FALSE, wrap.formula=" ~ key",
                       wrap.ncol=1, plot.type="box", x.lab="group",
                       y.lab=NULL, main=NULL,
                       compact.facet=TRUE, strip.position="right",
                       scales="free_y", hide.legend=TRUE, grid=FALSE,
                       text.size=14, text.x.angle=45, text.x.hjust=1,
                       ...) {
    plot.type <- match.arg(plot.type, c("box", "violin"))
    y.var <- as.character(as.list(as.formula(formula))[[2]])
    x.var <- as.character(as.list(as.formula(formula))[[3]])
    colour.var <- colour.var %||% x.var
    data[[colour.var]] <- factor(data[[colour.var]], levels = unique(data[[colour.var]]))
    p <- ggplot(data, aes_string(x = x.var, y = y.var, colour = colour.var))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (plot.type == "box") {
        p <- p + geom_boxplot()
    } else if (plot.type == "violin") {
        p <- p + geom_violin()
    }
    if (!is.null(x.lab)) p <- p + xlab(x.lab)
    if (!is.null(y.lab)) p <- p + ylab(y.lab)
    if (!is.null(main)) p <- p + ggtitle(main)
    if (compact.facet) {
        p <- p + theme(panel.spacing = unit(0, "lines"))
    }
    if (hide.legend) {
         p <- p + theme(legend.position = "none")
    }
    if (!grid) p <- p + theme(panel.grid = element_blank())
    p <- p + theme(text = element_text(size = text.size), axis.text.x = element_text(angle=text.x.angle, hjust=text.x.hjust))
    # Colouring setup
    nlevels <- length(levels(data[[colour.var]]))
    nc <- length(colour)
    n <-  floor(nlevels / nc)
    nmod <- nlevels %% nc
    if (nmod == 0) {
        p <- p + scale_colour_manual(values=rep(colour, n))
    } else {
        p <- p + scale_colour_manual(values=c(rep(colour, n), colour[1:n]))
    }
    p
}
##' @title boxplot.GRanges
##' @describeIn boxplot.pg Make boxplot of GRanges object
##' @param per.site normalize statistics to per-site values
##' @param ... arguments to pass to boxplot.pg
##' @export
boxplot.GRanges <- function(formula = "value ~ population", data=NULL, which=levels(factor(data$key)), per.site=TRUE, ...) {
    which <- match.arg(which, levels(factor(data$key)), several.ok = TRUE)
    df <- cbind(seqnames = as.character(seqnames(data)), as.data.frame(values(data)))
    df$ranges <- paste(as.character(seqnames(data)), start(data), "-", end(data), ":", sep = " ")
    y.var <- as.character(as.list(as.formula(formula))[[2]])
    if (per.site) df[[y.var]] <- df[[y.var]] / width(data)
    boxplot.pg(formula = formula, data = subset(df, key %in% which), ...)
}
##' @title boxplot.pg.detail
##' @describeIn boxplot.pg Make boxplot of detail statistiscs
##' @export
boxplot.pg.detail <- function(data, which=c("MDG1", "MDG2", "MDSD"), main="detail statistics", ...) {
    which = match.arg(which, c("MDG1", "MDG2", "MDSD"), several.ok = TRUE)
    boxplot.pg(data=subset(data, key %in% which), ...)
}
##' @title boxplot.pg.neutrality
##' @describeIn boxplot.pg Make boxplot of netrality statistics
##' @export
boxplot.pg.neutrality <- function(data=NULL, which=c("Tajima.D", "Fu.Li.F", "Fu.Li.F"), wrap=TRUE, ...) {
    which = match.arg(which, c("Tajima.D", "n.segregating.sites", "Rozas.R_2", "Fu.Li.F", "Fu.Li.D", "Fu.F_S", "Fay.Wu.H", "Zeng.E", "Strobeck.S"), several.ok = TRUE)
    boxplot.pg(data = subset(data, key %in% which), wrap = wrap, ...)
}
##' @title boxplot.pg.fixed
##' @describeIn boxplot.pg Make boxplot of fixed sites
##' @export
boxplot.pg.fixed <- function(data=NULL, main="Fixed sites", ...) {
    boxplot.pg(data = data, main = main, ...)
}
##' @title boxplot.pg.shared
##' @describeIn boxplot.pg Make boxplot of shared sites
##' @export
boxplot.pg.shared <- function(data=NULL, main="Shared sites", ...) {
    boxplot.pg(data = data, main = main, ...)
}
##' @title boxplot.pg.diversity
##' @describeIn boxplot.pg Make boxplot of diversity
##' @export
boxplot.pg.diversity <- function(data=NULL, which=c("nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"), wrap=TRUE, ...) {
    which <- match.arg(which, c("hap.diversity.within", "hap.F_ST.vs.all", "nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"), several.ok = TRUE)
    boxplot.pg(data = subset(data, key %in% which), wrap = wrap, ...)
}
##' @title boxplot.pg.diversity.between
##' @describeIn boxplot.pg Make boxplot of between population diversity
##' @export
boxplot.pg.diversity.between <- function(data=NULL, ...) {
    boxplot.pg(data = data, ...)
}

##' @title boxplot.pg.F_ST
##' @describeIn boxplot.pg Make boxplot of F_ST
##' @export
boxplot.pg.F_ST <- function(data=NULL, formula="value ~ seqnames", which=c("nucleotide.F_ST", "Nei.G_ST"), wrap=TRUE, ...) {
    which <- match.arg(which, c("haplotype.F_ST", "nucleotide.F_ST", "Nei.G_ST", "Hudson.G_ST", "Hudson.H_ST", "Hudson.K_ST"), several.ok = TRUE)
    boxplot.pg(formula = formula, data = subset(data, key %in% which), wrap = wrap, ...)
}
##' @title boxplot.pg.F_ST.pairwise
##' @describeIn boxplot.pg Make boxplot of pairwise F_ST
##' @export
boxplot.pg.F_ST.pairwise <- function(data=NULL, ...) {
    boxplot.pg(data = data, ...)
}

##' @title boxplot.pg.segregating.sites
##' @describeIn boxplot.pg Make boxplot of segregating sites
##' @export
boxplot.pg.segregating.sites <- function(data=NULL, ...) {
    boxplot.pg(data = data, ...)
}
