.application <- "Popoolation"
.varianceSlidingColnames <- c("seqnames", "position", "segregating.sites", "coverage", "score")
.varianceSlidingOutputColumns <- c("segregating.sites", "coverage", "score")

##' readVarianceSliding
##'
##' Read VarianceSliding results from popoolation
##'
##' @title
##'
##'
##' @return
##'
##' @author Per Unneberg
readVarianceSliding <- function(filename, measure="pi", ...) {
    message("read popoolation VarianceSliding output for file ", filename)
    measure <- match.arg(measure, c("pi", "D", "theta"))
    data <- read.table(filename)
    colnames(data) <- .varianceSlidingColnames
    data$score[data$score == "na"] <- NA
    data$score <- as.numeric(as.character(data$score))
    ManticoreDF(data, measurement.name = measure, application = .application)
}


##' VarianceSlidingAssays
##'
##' Read multiple VarianceSliding results from popoolation
##'
##' @param input.df Input data fram that must have columns filename, sample, and measure
##' @param colData column data that describes the samples
##' @param ... extra parameters passed to readVarianceSliding
##'
##' @return ManticoreRSE
##' @author Per Unneberg
##'
##' @importFrom GenomicRanges makeGRangesFromDataFrame
##' @importFrom SummarizedExperiment DataFrame
##' @importFrom BiocParallel bplapply
##'
VarianceSlidingAssays <- function(input.df, colData = NULL, width = NULL, ...) {
    stopifnot(inherits(input.df, c("data.frame", "DataFrame")))
    columns <- c("filename", "sample", "measure")
    stopifnot(columns %in% colnames(input.df))

    .loadData <- function(x) {
        measure <- unique(x$measure)
        data <- bplapply(as.character(x$filename),
                         function(y) {
            readVarianceSliding(y, measure = measure, ...)})
        names(data) <- x$sample
        data
    }
    rawData <- by(input.df, input.df$measure, .loadData)

    ## 1. create GRanges object as union of all sequence positions
    if (is.null(width))
        width <- rawData[[1]][[1]]$position[2] - rawData[[1]][[1]]$position[1]
    data <- data.frame(
        do.call("rbind", lapply(unlist(rawData),
                                  function(x) {
                             cbind(seqnames = as.character(x$seqnames), position = x$position)})))
    data$position <- as.numeric(as.character(data$position))
    data$start <- data$position - width / 2 + 1
    data$end <- data$position + width / 2
    gr <- unique(makeGRangesFromDataFrame(data))

    ## 2. modify raw data, adding NA elements to missing rows
    .makeGRanges <- function(y) {
        y$start <- y$position - width / 2 + 1
        y$end <- y$position + width / 2
        tmp <- makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
        i <- match(tmp, gr)
        x <- GRanges(gr)
        for (col in .varianceSlidingOutputColumns) {
            mcols(x)[[col]] <- rep(NA, length(gr))
            mcols(x)[[col]][i[!is.na(i)]] <- mcols(tmp)[[col]]
        }
        x
    }
    grlist <- lapply(rawData, function(x) {lapply(x, .makeGRanges)})
    .makeAssayData <- function(x) {
        ManticoreDF(do.call("cbind", lapply(grlist[[x]], function(y){mcols(y)})),
                    measurement.name = x,
                    application = .application)
    }
    assayData <- lapply(names(grlist), .makeAssayData)
    names(assayData) <- names(grlist)
    if (is.null(colData))
        colData = DataFrame(sample = input.df$sample, variable = .varianceSlidingOutputColumns)
    ManticoreRSE(assays = assayData,
                 rowRanges = gr, colData = colData)
}

