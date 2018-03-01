##' Get stats from a PopGenome GENOME object
##'
##' Retrieve genome stats for a PopGenome GENOME instance. Note that
##' this is only intended for use with some slots/statistical
##' operations. The PopGenome function does not work on data that has
##' been saved in rda files. They should be saved with save.session,
##' but this does not seem to work.
##'
##' @title getGenomeStats
##' @param object An R object
##' @param stats statistic to return
##' @param use.population.names use population names of object for plotting labels
##' @param use.region.names Use the region names as row names
##' @param long.format Return data frame in long output format
##' @param ... Arguments to pass to data access functions
##' @return A data frame in long or wide format, suitable for plotting with
##'     ggplot2
##' @author Per Unneberg
##' @export
setGeneric("getGenomeStats", function(object, stats=c("detail"), use.population.names=FALSE, use.region.names=FALSE, long.format=TRUE, ...) standardGeneric("getGenomeStats"))

##' @describeIn getGenomeStats Retrieve and concatenate data from a list of PopGenome GENOME instances.
##'
##' @examples
##' library(PopGenome)
##' # Read a vcf and generate two genome objects
##' vcf_file <- system.file("extdata", "medium.call.vcf.gz")
##' scaffold1 <- readVCF(vcf_file, 1000, frompos=1, topos=1000000, tid="scaffold1")
##' scaffold13 <- readVCF(vcf_file, 1000, frompos=1, topos=1000000, tid="scaffold13")
##' slist <- list(scaffold1=scaffold1, scaffold13=scaffold13)
##'
setMethod("getGenomeStats", "list", function(object, stats, use.population.names, use.region.names, long.format, ...) {
    if (!requireNamespace("PopGenome", quietly = TRUE)) {
        stop("Package \"PopGenome\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("Package \"tidyr\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    stopifnot(all(unlist(lapply(object, function(x){inherits(x, "GENOME")}))))
    stats <- match.arg(stats, c("summary", "detail", "neutrality",
                                "fixed.shared", "diversity",
                                "diversity.between", "F_ST",
                                "F_ST.pairwise",
                                "linkage", "sweeps", "recomb"))
    res <- data.frame()
    i <- 1
    populations <- names(object[[1]]@populations)
    pairs <- combn(populations, 2, function(x){paste(x, collapse="/")})
    fun <- list(summary = get.sum.data, detail = get.detail, neutrality = get.neutrality,
                fixed.shared = get.fixed.shared, diversity =
                get.diversity, diversity.between =
                get.diversity.between, F_ST = get.F_ST, F_ST.pairwise
                = get.F_ST.pairwise, linkage = get.linkage, sweeps =
                get.sweeps, recomb = get.recomb)

    f <- fun[[stats]]
    gather.key = "key"
    for (name in names(object)) {
        regions <- object[[name]]@region.names
        if (stats %in% c("detail", "neutrality", "diversity", "linkage", "sweeps", "recomb")) {
            tmp <- f(object[[name]])
            if (use.population.names) {
                rownames(tmp) <- populations
            }
            pos <- seq(i, i + length(regions) - 1)
            tmp <- do.call("rbind", lapply(rownames(tmp), function(x){data.frame(population=x, region=rownames(tmp[x,][[1]]), tmp[x,][[1]], pos = pos)}))
            tmp$name <- name
            rownames(tmp) <- NULL
            gather.exclude <- c("population", "region", "name", "pos")
        } else if (stats %in% c("fixed.shared", "diversity.between")) {
            tmp <- as.data.frame(f(object[[name]], ...))
            if (use.population.names) {
                colnames(tmp) <- pairs
            }
            tmp$name <- name
            tmp$region <- regions
            tmp$pos <- seq(i, i + length(regions) - 1)
            gather.key <- "population"
            gather.exclude <- c("region", "name", "pos")
        } else if (stats %in% c("F_ST.pairwise")) {
            tmp <- f(object[[name]])
            col.names <- colnames(tmp[[1]])
            pos <- seq(i, i + length(regions) - 1)
            tmp <- do.call(
                "rbind", lapply(rownames(tmp),
                                function(x) {
                             data.frame(key = x, region = regions, pos = pos, tmp[x, ])}))
            if (use.population.names) {
                colnames(tmp)[4:dim(tmp)[2]] <- pairs
            } else {
                colnames(tmp)[4:dim(tmp)[2]] <- col.names
            }
            tmp$name <- name
            gather.key <- "population"
            gather.exclude <- c("key", "region", "pos", "name")
        } else if (stats %in% c("summary")) {
            message("Analyzing ", stats, " data")
            ## No population data here; just summary
            ## Result is a matrix that can readily be converted to a data frame
            tmp <- as.data.frame(f(object[[name]]))
            tmp$name <- name
            gather.exclude <- c("name")
        } else if (stats %in% c("F_ST")) {
            tmp <- as.data.frame(f(object[[name]]))
            tmp$region <- regions
            tmp$name <- name
            tmp$pos <- seq(i, i + length(regions) - 1)
            gather.exclude <- c("region", "name", "pos")
        } else {
            stop("shouldn't end up here")
        }
        i <- i + length(regions)
        message("Adding data for ", name)
        res <- rbind(res, tmp)
    }
    if (long.format) {
        res <- gather(res, gather.key, "value", -gather.exclude)
        colnames(res)[which(colnames(res) == "gather.key")] <- gather.key
    }
    class(res) <- c(paste("pg", stats, sep="."), "data.frame")
    return(res)
})



## Write new functions for each case, and a function for annotating
## data, as in bioodo
setGeneric("get.fixed.shared", function(object, which="fixed", ...) standardGeneric("get.fixed.shared"))
setMethod("get.fixed.shared", "GENOME", function(object, which, ...) {
    which <- match.arg(which, c("fixed", "shared"))
    message("Getting ", which, " counts")
    slots <- list(fixed = "n.fixed.sites", shared = "n.shared.sites")
    return (slot(object, slots[[which]]))
})

setGeneric("get.diversity.between", function(object, which="nuc", ...) standardGeneric("get.diversity.between"))
setMethod("get.diversity.between", "GENOME", function(object, which, ...) {
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
})

setGeneric("get.F_ST.pairwise", function(object, ...) standardGeneric("get.F_ST.pairwise"))
setMethod("get.F_ST.pairwise", "GENOME", function(object, which, ...) {
    df <- get.F_ST(object, pairwise = TRUE)
    return (df)
})

##' genomewide.stats
##'
##' Calculate genome wide statistics
##' @title genomewide.stats
##' @param object GENOME object on which to perform calculations
##' @param which statistics to calculate
##' @param biallelic.structure calculate biallelic structure
##' @param ...
##' @return Updated object
##' @author Per Unneberg
##' @export
setGeneric("genomewide.stats", function(object, which=c("detail", "neutrality", "F_ST", "diversity", "diversity.between", "fixed.shared", "linkage"), biallelic.structure=TRUE, ...) standardGeneric("genomewide.stats"))

setMethod("genomewide.stats", "GENOME", function(object, which, ...) {
    which <- match.arg(which, c("detail", "neutrality", "F_ST", "diversity", "diversity.between", "fixed.shared", "linkage", "R2", "recomb", "sweeps"), several.ok=TRUE)
    if ("detail" %in% which) {
        object <- detail.stats(object, biallelic.structure=biallelic.structure, ...)
    }
    if ("neutrality" %in% which) {
        object <- neutrality.stats(object, ...)
    }
    if ("F_ST" %in% which) {
        object <- F_ST.stats(object, ...)
    }
    if ("diversity" %in% which) {
        object <- diversity.stats(object, ...)
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


##' plot.pg
##'
##' Generic plotting function for PopGenome results
##' @title plot.pg
##' @param data long format from getGenomeStats function
##' @param ...
##' @return ggplot
##' @author Per Unneberg
plot.pg <- function(data, x="pos", y="value", colour="name", wrap=TRUE, wrap.formula="key ~ population", wrap.ncol=1, plot.type="point",
                    x.lab=NULL, y.lab=NULL, main=NULL, colour.scale=NULL, compact.facet=TRUE,
                    strip.position="right", scales="free_y", ...) {
    plot.type <- match.arg(plot.type, c("point", "line"))
    p <- ggplot(data, aes_string(x = x, y = y, colour = colour))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (plot.type == "point") {
        p <- p + geom_point()
    } else if (plot.type == "line") {
        p <- p + geom_line()
    }
    if (!is.null(x.lab)) p <- p + xlab(x.lab)
    if (!is.null(y.lab)) p <- p + ylab(y.lab)
    if (!is.null(main)) p <- p + ggtitle(main)
    if (compact.facet) {
        p <- p + theme(panel.spacing = unit(0, "lines"))
    }
    p
}

##' @export
plot.pg.summary <- function(data, x="name", y="value", wrap.formula="~ key", which=c("n.sites", "trans.transv.ratio"), ...) {
    which <- match.arg(which, c("n.sites", "n.biallelic.sites", "n.gaps", "n.unknowns", "n.valid.sites", "n.polyallelic.sites", "trans.transv.ratio"), several.ok=TRUE)
    plot.pg(subset(data, key %in% which), x = x, y = y, wrap.formula = wrap.formula, ...)
}

##' @export
plot.pg.detail <- function(data, which=c("MDG1", "MDG2", "MDSD"), ...) {
    which = match.arg(which, c("MDG1", "MDG2", "MDSD"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}

##' @export
plot.pg.neutrality <- function(data, which=c("Tajima.D", "Fu.Li.F", "Fu.Li.F"), ...) {
    which = match.arg(which, c("Tajima.D", "n.segregating.sites", "Rozas.R_2", "Fu.Li.F", "Fu.Li.D", "Fu.F_S", "Fay.Wu.H", "Zeng.E", "Strobeck.S"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}

##' @export
plot.pg.fixed.shared <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}

##' @export
plot.pg.diversity <- function(data, which=c("nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"),  ...) {
    which <- match.arg(which, c("hap.diversity.within", "hap.F_ST.vs.all", "nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"), several.ok=TRUE)
    plot.pg(subset(data, key %in% which), ...)
}

##' @export
plot.pg.diversity.between <- function(data, wrap.formula="~ population", colour="name", ...) {
    plot.pg(data, wrap.formula = wrap.formula, colour = colour, ...)
}

##' @export
plot.pg.F_ST <- function(data, wrap.formula="~ key", which=c("nucleotide.F_ST", "Nei.G_ST"), ...) {
    which <- match.arg(which, c("haplotype.F_ST", "nucleotide.F_ST", "Nei.G_ST", "Hudson.G_ST", "Hudson.H_ST", "Hudson.K_ST"), several.ok=TRUE)
    plot.pg(subset(data, key %in% which), wrap.formula = wrap.formula, ...)
}

##' @export
plot.pg.F_ST.pairwise <- function(data, ...) {
    plot.pg(data, ...)
}

