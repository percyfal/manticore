##' BcftoolsStats object
##'
##' Representation of bcftools stats output
##'
##' @title BcftoolsStats
##' @slot ID ID section
##' @slot SN summary numbers section
##' @slot TSTV transition/transversion section
##' @slot SiS singleton stats
##' @slot AF stats by non-reference allele frequency
##' @slot QUAL stats by quality
##' @slot IDD indel distribution
##' @slot ST substitution types
##' @slot DP depth distribution
##'
##' @importFrom methods callNextMethod
##'
BcftoolsStats<- setClass("BcftoolsStats", slots=c(ID="data.frame", SN="data.frame",
                                                  TSTV="data.frame",
                                                  SiS="data.frame", AF="data.frame",
                                                  QUAL="data.frame",
                                                  IDD="data.frame", ST="data.frame",
                                                  DP="data.frame", label="character"))

setMethod("initialize", "BcftoolsStats", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})

##' readBcftoolsStats
##'
##' Read output from bcftools stats
##'
##' @title readBcftoolsStats
##' @param filename  File to read
##' @param label Label to assign object (e.g. sample)
##' @return BcftoolsStats object
##' @author Per Unneberg
##'
##' @export
##'
##' @importFrom methods slot slot<- slotNames
##' @importFrom utils read.delim
##'
setGeneric("readBcftoolsStats", function(filename, label=NULL) {standardGeneric("readBcftoolsStats")})

setMethod("readBcftoolsStats", signature="character", definition=function(filename, label) {
    con <- file(filename, open = "r")
    lines <- readLines(con)
    # Loop the labels and convert to data frames
    obj <- BcftoolsStats()
    for (x in slotNames(obj)) {
        if (x == "label") next
        header <- gsub(
            "\\[[0-9]+\\]", "",
            unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep = ""), lines))]),
                            "\t")))
        slot(obj, x, check = TRUE) <- read.delim(
            tc <- textConnection(lines[grepl(paste(x, "\\t", sep = ""), lines)]),
            header = FALSE, skip = 1, col.names = header)
        close(tc)
        slot(obj, x)[x] <- NULL
    }
    close(con)
    if (is.null(label)) {
        obj@label <- filename
    } else {
        obj@label <- label
    }
    obj
})


##' read.bcftools.stats
##'
##' Read output from bcftools stats
##'
##' @title read.bcftools.stats
##' @param filename File to read
##' @param label Label to assign object (e.g. sample)
##' @return bcftools.stats object
##' @author Per Unneberg
##' @export
read.bcftools.stats <- function(filename, label=NULL) {
    con <- file(filename, open = "r")
    lines <- readLines(con)
    sections <- c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", "label")
    obj <- list()
    # Loop the labels and convert to data frames
    for (x in sections) {
        if (x == "label") next
        header <- gsub("\\[[0-9]+\\]", "", unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep=""), lines))]), "\t")))
        obj[[x]] <- read.delim(tc <- textConnection(lines[grepl(paste(x, "\\t", sep = ""), lines)]),
                               header = FALSE, skip = 1, col.names = header)
        close(tc)
        obj[[x]][x] <- NULL
    }
    close(con)
    if (is.null(label)) {
        message("setting label to filename")
        obj[["label"]] <- filename
    } else {
        obj[["label"]] <- label
    }
    class(obj) <- "bcftools.stats"
    obj
}

##' summary.bcftools.stats
##'
##' Print summary for bcftools stats
##'
##' @title summary.bcftools.stats
##' @param obj bcftools.stats object
##' @return data frame summary
##' @author Per Unneberg
##' @export
summary.bcftools.stats <- function(obj) {
    as.data.frame(obj$SN)
}


##' gplot.bcftools.stats
##'
##' Plot summary of bcftools stats results
##' @title plot.bcftools.stats
##' @param obj bcftools.stats object
##' @param which which plots to produce (defalt all)
##' @param ncol number of columns in grid plot
##' @param text.size text size of plots
##' @param theme.default default theme (theme_bw)
##' @param ... parameters passed to generic plot function
##' @return ggplot2 object, or a plot
##' @author Per Unneberg
##' @export
##'
##' @import tidyr
##'
gplot.bcftools.stats <- function(obj, which=c("SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP"), ncol=2, text.size=12, theme.default=theme_bw(), ...) {
    which <- match.arg(which, c("SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP"), several.ok=TRUE)
    message("Producing ", length(which), " plots")
    old <- theme_set(theme.default)
    theme_update(text = element_text(size = text.size))
    plist <- list()
    if ("SN" %in% which) {
        data <- obj$SN
        plist$SN <- ggplot(data = data, aes(y = key, x = value)) + geom_point(...) + ggtitle("Summary numbers") + theme(text = element_text(size = text.size))
    }
    if ("TSTV" %in% which) {
        data <- gather(obj$TSTV, key, value, -id)
        tstv <- c("ts.tv", "ts.tv..1st.ALT.")
        plist$TSTV <- arrangeGrob(
            ggplot(data = subset(data, key %in% tstv), aes(y = key, x = value)) + geom_point(...) + ggtitle("Ti/Tv") + theme(text = element_text(size = text.size)),
            ggplot(data = subset(data, !(key %in% tstv)), aes(y = key, x = value)) + geom_point(...) + ggtitle("Ti/Tv, counts") + theme(text = element_text(size = text.size)),
            nrow = 2)
    }
    if ("SiS" %in% which) {
        data <- gather(obj$SiS, key, value, -id)
        plist$SiS <- ggplot(data = data, aes(y = key, x = value)) + geom_point(...) + ggtitle("Singleton stats") + theme(text = element_text(size = text.size))
    }
    if ("AF" %in% which) {
        data <- gather(obj$AF, key, value, -allele.frequency, -id)
        plist$AF <- ggplot(data = data, aes(x = allele.frequency, y = value, color = key)) + geom_point(...) + ggtitle("Stats non-reference allele frequency") + theme(text = element_text(size = text.size))
    }
    if ("QUAL" %in% which) {
        data <- gather(obj$QUAL, key, value, -Quality, -id)
        plist$QUAL <- ggplot(data = data, aes(x = Quality, y = value, color = key)) + geom_point(...) + ggtitle("Stats by quality") + theme(text = element_text(size = text.size))
    }
    if ("IDD" %in% which) {
        data <- obj$IDD
        plist$IDD <- ggplot(data = data, aes(x = length..deletions.negative., y = count)) + geom_point(...) + ggtitle("InDel distribution") + xlab("Length (deletions negative)") + theme(text = element_text(size = text.size))
    }
    if ("ST" %in% which) {
        data <- obj$ST
        plist$ST <- ggplot(data = data, aes(x = type, y = count)) + geom_point(...) + ggtitle("Substitution types") + theme(text = element_text(size = text.size))
    }
    if ("DP" %in% which) {
        data <- gather(obj$DP, key, value, -bin, -id)
        plist$DP <- ggplot(data = data, aes(x = bin, y = value, color = key)) + geom_point(...) + ggtitle("Depth distribution") + theme(text = element_text(size = text.size))
    }

    if (length(plist) > 1) {
        p <- arrangeGrob(grobs = plist, ncol = ncol)
    } else {
        p <- plist[[1]]
    }
    theme_set(old)
    p
}
