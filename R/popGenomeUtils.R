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
                         "outgroup", "region.names", "feature.names",
                         "genelength", "n.sites", "n.sites2",
                         "n.biallelic.sites", "n.gaps", "n.unknowns",
                         "n.valid.sites", "n.polyallelic.sites",
                         "trans.transv.ratio", "keep.start.pos",
                                                  "Coding.region", "UTR.region",
                         "Intron.region", "Exon.region",
                         "Gene.region", "Pop_Neutrality", "Pop_FSTN",
                         "Pop_FSTH", "Pop_Linkage", "Pop_Slide",
                         "Pop_MK", "Pop_Detail", "Pop_Recomb",
                         "Pop_Sweeps", "FSTNLISTE",
                             "nucleotide.F_ST2",
                                   "region.data", "region.stats")
.region.stats.slots <- c()

##' Get data from a PopGenome GENOME object
##'
##' Retrieve data for a specific slot in a PopGenome GENOME instance.
##' Note that this is only intended for use with some slots.
##'
##' @title getGenomeData
##' @param object An R object
##' @param slotName The name of the slot to retrieve
##' @param population.names Rename population columns in GENOME
##' @param use.region.names Use the region names as row names
##' @return A data frame in long format, suitable for plotting with ggplot2
##' @author Per Unneberg
##' @export
##'
##'

setGeneric("test", function(tmp) standardGeneric("test"))


##setGeneric("getGenomeData", function(object) 1 ) #, function(object, slotName, col.names=NULL, use.region.names=FALSE, ...) standardGeneric("getGenomeData"))


##' @describeIn getGenomeData Retrieve and concatenate data from a list of PopGenome GENOME instances.
##' @import PopGenome
##'
##' @examples
##' library(PopGenome)
##' # Read a vcf and generate two genome objects
##' scaffold_1 <- readVCF(vcf_file, 1000, frompos=1, topos=1000, tid="scaffold_1")
##' scaffold_2 <- readVCF(vcf_file, 1000, frompos=1, topos=1000, tid="scaffold_2")
##' slist <- list(scaffold_1=scaffold_1, scaffold_2=scaffold_2)
##' df.segregating.sites <- getGenomeData(slist, "n.segregating.sites"

## setMethod("getGenomeData", "list", function(object, slotName, load.data, col.names, use.region.names, ...) {
##     if (!requireNamespace("PopGenome", quietly = TRUE)) {
##         stop("Package \"PopGenome\" needed for this function to work. Please install it.",
##              call. = FALSE)
##     }
##     stopifnot(all(unlist(lapply(slist, function(x){inherits(x, "GENOME")}))))
##     df = data.frame()
##     i = 1
##     for (oname in names(object)) {
##         tmp <- as.data.frame(slot(object[[oname]], slotName))
##         if (!is.null(col.names)) {
##             colnames(tmp) <- col.names
##         }
##         df2 <- cbind(gather(tmp), oname)
##         if (use.region.names) {
##             rownames(df2) <- slot(object[[oname]], "region.names")
##         }
##         df2$pos <- rep(seq(i, i + dim(tmp)[1] - 1), dim(tmp)[2])
##         message("Adding data for ", oname)
##         df <- rbind(df, df2)
##         i <- i + dim(tmp)[1]
##     }
##     df
## })
