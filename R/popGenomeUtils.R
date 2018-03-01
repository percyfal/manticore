.GENOME.class.slots <- c(
    "nucleotide.F_ST",                  # get.F_ST, get.diversity, F_ST.stats; get.F_ST gets overall stats, get.diversity between-population stats
    "nuc.diversity.between",             #
    "nuc.diversity.within",
    "nuc.F_ST.pairwise",
    "nuc.F_ST.vs.all",
    "n.haplotypes",
    "hap.diversity.within",
    "hap.diversity.between",
    "Pi",
    "PIA_nei",
    "haplotype.counts",
    "haplotype.F_ST",
    "hap.F_ST.pairwise",
    "Nei.G_ST.pairwise",
    "hap.F_ST.vs.all",
    "Nei.G_ST",
    "Hudson.G_ST",
    "Hudson.H_ST",
    "Hudson.K_ST",
    "Hudson.Snn",
    "Phi_ST",
    "hap.pair.F_ST",
    "MKT",
    "Tajima.D",
    "SLIDE",
    "Fay.Wu.H",
    "Zeng.E",
    "theta_Tajima",
    "theta_Watterson",
    "theta_Fu.Li",
    "theta_Achaz.Watterson",
    "theta_Achaz.Tajima",
    "theta_Fay.Wu",
    "theta_Zeng",
    "Fu.Li.F",
    "Fu.Li.D",
    "Yach",
    "n.segregating.sites",
    "Rozas.R_2",
    "Fu.F_S",
    "Strobeck.S",
    "Kelly.Z_nS",
    "Rozas.ZZ",
    "Rozas.ZA",
    "Wall.B",
    "Wall.Q",
    "mult.Linkage",
    "RM",
    "CL",
    "CLmax",
    "CLR",
    "MDSD",
    "MDG1",
    "MDG2",
    "D",
    "BD",
    "BDF",
    "BDF_bayes",
    "alpha_ABBA",
    "alpha_BABA",
    "beta_BBAA",
    "Bd_clr",
    "Bd_dir",
    "P.Bd_clr",
    "f",
    "RNDmin",
    "D.z",
    "D.pval",
    "jack.knife",
    "missing.freqs",
    "n.fixed.sites",
    "n.shared.sites",
    "n.monomorphic.sites",
    "genes"
)
.GENOME.class.slots.blacklist <- c("BIG.BIAL", "SLIDE.POS",
                                   "big.data","gff.info", "snp.data",
                                   "basepath", "project",
                                   "populations", "poppairs",
                                   "outgroup", "region.names",
                                   "feature.names", "genelength",
                                   "n.sites", "n.sites2",
                                   "n.biallelic.sites", "n.gaps",
                                   "n.unknowns", "n.valid.sites",
                                   "n.polyallelic.sites",
                                   "trans.transv.ratio",
                                   "keep.start.pos", "Coding.region",
                                   "UTR.region", "Intron.region",
                                   "Exon.region", "Gene.region",
                                   "Pop_Neutrality", "Pop_FSTN",
                                   "Pop_FSTH", "Pop_Linkage",
                                   "Pop_Slide", "Pop_MK",
                                   "Pop_Detail", "Pop_Recomb",
                                   "Pop_Sweeps", "FSTNLISTE",
                                   "nucleotide.F_ST2", "region.data",
                                   "region.stats")
.region.stats.slots <- c()

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
##' @param ...
##' @return A data frame in long format, suitable for plotting with
##'     ggplot2
##' @author Per Unneberg
##' @export
setGeneric("getGenomeStats", function(object, stats=c("detail"), use.population.names=FALSE, use.region.names=FALSE, ...) standardGeneric("getGenomeStats"))

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
setMethod("getGenomeStats", "list", function(object, stats, use.population.names, use.region.names, ...) {
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
    for (name in names(object)) {
        regions <- object[[name]]@region.names
        if (stats %in% c("detail", "neutrality", "diversity", "linkage", "sweeps", "recomb")) {
            tmp <- f(object[[name]])
            if (use.population.names) {
                rownames(tmp) <- populations
            }
            tmp <- do.call("rbind", lapply(rownames(tmp), function(x){data.frame(population=x, region=rownames(tmp[x,][[1]]), tmp[x,][[1]])}))
            message("adding plot positions")
            tmp$pos <- rep(seq(i, i + length(regions) - 1), length(populations))
            i <- i + length(regions)
            tmp$name <- name
            tmp <- gather(tmp, key, value, -population, -region, -name, -pos)
        } else if (stats %in% c("fixed.shared", "diversity.between")) {
            tmp <- as.data.frame(f(object[[name]], ...))
            if (use.population.names) {
                colnames(tmp) <- pairs
            }
            tmp <- gather(tmp, population, value)
            tmp$name <- name
            tmp$region <- regions
            tmp$pos <- seq(i, i + length(regions) - 1)
            i <- i + length(regions)
        } else if (stats %in% c("F_ST.pairwise")) {
            tmp <- f(object[[name]])
            col.names <- colnames(tmp[[1]])
            pos <- seq(i, i + length(regions) - 1)
            i <- i + length(regions)
            tmp <- do.call(
                "rbind", lapply(rownames(tmp),
                                function(x) {
                             data.frame(key = x, region = regions, pos = pos, tmp[x, ])}))
            if (use.population.names) {
                colnames(tmp)[4:dim(tmp)[2]] <- pairs
            } else {
                colnames(tmp)[4:dim(tmp)[2]] <- col.names
            }
            tmp <- gather(tmp, population, value, -key, -region, -pos)
            tmp$name <- name
        } else if (stats %in% c("summary")) {
            message("Analyzing ", stats, " data")
            ## No population data here; just summary
            ## Result is a matrix that can readily be converted to a data frame
            tmp <- as.data.frame(f(object[[name]]))
        } else if (stats %in% c("F_ST")) {
            tmp <- as.data.frame(f(object[[name]]))
            tmp <- gather(tmp)
            tmp$region <- regions
            tmp$name <- name
            tmp$pos <- seq(i, i + length(regions) - 1)
            i <- i + length(regions)
        } else {
            stop("shouldn't end up here")
        }
        message("Adding data for ", name)
        res <- rbind(res, tmp)
    }
    class(res) <- c(stats, "data.frame")
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
