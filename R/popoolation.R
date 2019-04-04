.application <- "Popoolation"
.varianceSlidingColnames <- c("seqnames", "position", "segregating.sites", "coverage", "score")
.varianceSlidingOutputColumns <- c("segregating.sites", "coverage", "score")


.readVarianceSlidingRaw <- function(filename, measure) {
    data <- read.table(filename)
    colnames(data) <- .varianceSlidingColnames
    data$score[data$score == "na"] <- NA
    data$score <- as.numeric(as.character(data$score))
    data
}


## Convert data.frame to ManticoreRSE object
.asManticoreRSE <- function(data, window.size, sample, measure) {
    if (is.null(window.size)) {
        window.size <- data$position[2] - data$position[1]
        message("window.size parameter undefined; inferring window size to ", window.size, " from data")
    }
    assayData <- list(DataFrame(score=data$score))
    names(assayData) <- measure
    data$start <- data$position - window.size / 2 + 1
    data$end <- data$position + window.size / 2
    sw <- SWindows(seqnames = data$seqnames, ranges = IRanges(start = data$start, end = data$end),
                   coverage = data$coverage, segregating.sites = data$segregating.sites,
                   window.size = window.size)
    colData = S4Vectors::DataFrame(sample = sample)
    ManticoreRSE(assays = assayData,
                 rowRanges = sw, colData = colData, window.size = window.size)
}

##' readVarianceSliding
##'
##' Read VarianceSliding results from popoolation. This function reads
##' data, adds start and end columns and converts data to a
##' ManticoreRSE object.
##'
##'
##' @param filename filename to parse
##' @param window.size window size
##' @param measure measure
##'
##' @return
##'
##' @author Per Unneberg
##'
readVarianceSliding <- function(filename, measure = "pi", window.size = NULL, sample) {
    measure <- match.arg(measure, c("pi", "D", "theta"))
    data <- .readVarianceSlidingRaw(filename, measure)
    .asManticoreRSE(data, window.size, sample, measure)
}


##' VarianceSlidingAssays
##'
##' Read multiple VarianceSliding results from popoolation
##'
##' @param input.df Input data fram that must have columns filename,
##'     sample, and measure
##' @param colData column data that describes the samples
##' @param window.size window size; if NULL, entire chromosomes are
##'     implied
##' @param class arbitrary classification of windows, e.g. to
##'     autosomes, sex chromosomes etc.
##' @param ... extra parameters passed to readVarianceSliding
##'
##' @return ManticoreRSE Assays consist of the requested statistics
##' @author Per Unneberg
##'
##' @importFrom GenomicRanges makeGRangesFromDataFrame, GenomicRangesList
##' @importFrom S4Vectors DataFrame
##' @importFrom BiocParallel bplapply
##'
VarianceSlidingAssays <- function(input.df, colData = NULL, window.size = NULL, class = NA, ...) {
    stopifnot(inherits(input.df, c("data.frame", "DataFrame")))
    columns <- c("filename", "sample", "measure")
    stopifnot(columns %in% colnames(input.df))

    .loadData <- function(x) {
        l <- as.list(by(x, x$measure, list))
        data <- BiocParallel::bplapply(l,
                                       function(y) {
                                  readVarianceSliding(as.character(y$filename),
                                                      measure = as.character(y$measure),
                                                      window.size,
                                                      sample = as.character(y$sample),
                                                      ...)})
        obj <- data[[1]]
        for (y in names(data))
            assay(obj, y) <- assay(data[[y]])
        colnames(obj) <- "score"
        obj
    }
    rawData <- lapply(as.list(by(input.df, input.df$sample, list)), .loadData)
    do.call(cbind, rawData)
}

