.application <- "Popoolation"
.varianceSlidingColnames <- c("seqnames", "position", "segregating.sites", "coverage", "score")
.varianceSlidingOutputColumns <- c("segregating.sites", "coverage", "score")


.readVarianceSlidingRaw <- function(filename) {
    data <- read.table(filename)
    colnames(data) <- .varianceSlidingColnames
    data$score[data$score == "na"] <- NA
    data$score <- as.numeric(as.character(data$score))
    data
}


## Convert data.frame to WindowedSummarizedExperiment object
.asWindowedSummarizedExperiment <- function(data, window.size, sample, measure) {
    if (is.null(window.size)) {
        window.size <- data$position[2] - data$position[1]
        message("window.size parameter undefined; inferring window size to ", window.size, " from data")
    }
    df <- DataFrame(score=data$score)
    colnames(df) <- sample
    assayData <- S4Vectors::SimpleList(df)
    names(assayData) <- measure
    data$start <- data$position - window.size / 2 + 1
    data$end <- data$position + window.size / 2
    sw <- SWindows(seqnames = data$seqnames, ranges = IRanges::IRanges(start = data$start, end = data$end),
                   coverage = data$coverage, segregating.sites = data$segregating.sites,
                   window.size = window.size)
    colData <- S4Vectors::DataFrame(sample = sample)
    WindowedSummarizedExperiment(assays = assayData,
                 rowRanges = sw, colData = colData)
}

##' readVarianceSliding
##'
##' Read VarianceSliding results from popoolation. This function reads
##' data, adds start and end columns and converts data to a
##' WindowedSummarizedExperiment object.
##'
##'
##' @param filename filename to parse
##' @param window.size window size
##' @param measure measure
##'
##' @examples
##' fn <- system.file("extdata", "popoolation", "dmel.A.D.txt.gz", package = "manticore")
##' readVarianceSliding(fn, measure = "D", sample = "A")
##'
##' @return WindowedSummarizedExperiment object of one assay
##'
##' @author Per Unneberg
##'
readVarianceSliding <- function(filename, measure = "pi", window.size = NULL, sample) {
    measure <- match.arg(measure, c("pi", "D", "theta"))
    data <- .readVarianceSlidingRaw(filename)
    .asWindowedSummarizedExperiment(data, window.size, sample, measure)
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
##' @return WindowedSummarizedExperiment Assays consist of the requested statistics
##' @author Per Unneberg
##'
##' @importFrom S4Vectors DataFrame
##' @importFrom BiocParallel bplapply
##'
##' @examples
##' dmel.df <- data.frame(filename = list.files(system.file("extdata", "popoolation", package = "manticore"), full.names = TRUE),
##'                       sample = c(rep("A", 3), rep("B", 3)), measure = rep(c("D", "pi", "theta"), 2))
##' mrse <- VarianceSlidingAssays(dmel.df)
##'
##'
VarianceSlidingAssays <- function(input.df, colData = NULL, window.size = NULL, class = NA, ...) {
    stopifnot(inherits(input.df, c("data.frame", "DataFrame")))
    columns <- c("filename", "sample", "measure")
    stopifnot(columns %in% colnames(input.df))
    input.df <- droplevels(input.df)

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
        obj
    }
    rawData <- lapply(as.list(by(input.df, input.df$sample, list)), .loadData)
    do.call(cbind, rawData)
}

