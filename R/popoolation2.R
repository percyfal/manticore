.application <- "Popoolation2"


.readPopoolation2.fet.fst <- function(filename, assay, window.size = NULL) {
    data <- read.table(filename)
    df <- as.data.frame(sapply(data[6:ncol(data)], as.character), stringsAsFactors = FALSE)
    data <- data[1:(ncol(data) - 1)]
    colnames(data) <- c("seqnames", "position", "segregating.sites", "coverage", "avmincov")
    if (is.null(window.size)) {
        window.size <- data$position[2] - data$position[1]
        message("window.size parameter undefined; inferring window size to ", window.size, " from data")
    }
    sample <- gsub("=.+", "", df[1,])
    df <- DataFrame(sapply(df, function(x){as.numeric(gsub(".+=", "", x))}))
    colnames(df) <- sample
    assayData <- list(df)
    names(assayData) <- assay
    data$start <- data$position - window.size / 2 + 1
    data$end <- data$position + window.size / 2
    sw <- SWindows(seqnames = data$seqnames, ranges = IRanges(start = data$start, end = data$end),
                   coverage = data$coverage, segregating.sites = data$segregating.sites,
                   window.size = window.size)
    colData <- S4Vectors::DataFrame(sample = sample)
    WindowedSummarizedExperiment(assays = assayData,
                                 rowRanges = sw, colData = colData)
}

.readPopoolation2.pwc.rc <- function(filename, assay) {
    if (assay == "pwc")
        n.mcol <- 8
    else
        n.mcol <- 9
    header <- unlist(strsplit(gsub("##", "", readLines(filename, n=1)), "\t"))
    data <- read.table(filename)
    colnames(data) <- header
    colnames(data)[1:2] <- c("seqnames", "position")
    df <- DataFrame(data[(n.mcol + 1):ncol(data)])
    colnames(df) <- header[(n.mcol + 1):length(header)]
    assayData <- list(df)
    names(assayData) <- assay
    gr <- GRanges(seqnames = data$seqnames, ranges = IRanges(start = data$position, end = data$position))
    mcols(gr) <- data[3:n.mcol]
    colData <- S4Vectors::DataFrame(colnames = header[(n.mcol + 1):length(header)])
    SummarizedExperiment(assays = assayData,
                         rowRanges = gr,
                         colData = colData)
}

##' @rdname readPopoolation2
##'
##' @title readPopoolation2
##'
##' @description Read popoolation2 results
##'
##' @param filename filename
##' @param assay assay type
##' @param window.size window size
##'
##' @return A WindowedSummarizedExperiment or RangedSummarizedExperiment, depending on assay
##' @author Per Unneberg
##'
##'
readPopoolation2 <- function(filename, assay = "fet", window.size = NULL) {
    assay <- match.arg(assay, c("fet", "fst", "cmh", "pwc", "rc"))
    .parsers <- list(fet = .readPopoolation2.fet.fst,
                     fst = .readPopoolation2.fet.fst,
                     pwc = .readPopoolation2.pwc.rc,
                     rc = .readPopoolation2.pwc.rc)
    if (assay %in% c("fet", "fst")) {
        return(.readPopoolation2.fet.fst(filename, assay, window.size))
    } else if (assay %in% c("rc", "pwc")) {
        return(.readPopoolation2.pwc.rc(filename, assay))
    }
}
