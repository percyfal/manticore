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
setGeneric("readBcftoolsStats", function(filename, label=NULL) {standardGeneric("readBcftoolsStats")})

setMethod("readBcftoolsStats", signature="character", definition=function(filename, label) {
    con <- file(filename, open = "r")
    lines <- readLines(con)
    # Loop the labels and convert to data frames
    obj <- BcftoolsStats()
    for (x in slotNames(obj)) {
        if (x == "label") next
        header <- gsub("\\[[0-9]+\\]", "", unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep=""), lines))]), "\t")))
        slot(obj, x, check = TRUE) <- read.delim(textConnection(lines[grepl(paste(x, "\\t", sep=""), lines)]), header=FALSE, skip=1, col.names=header)
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
read.bcftools.stats <- function(filename, label=NULL) {
    con <- file(filename, open = "r")
    lines <- readLines(con)
    sections <- c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", "label")
    obj <- list()
    # Loop the labels and convert to data frames
    for (x in sections) {
        if (x=="label") next
        header <- gsub("\\[[0-9]+\\]", "", unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep=""), lines))]), "\t")))
        obj[[x]] <- read.delim(textConnection(lines[grepl(paste(x, "\\t", sep=""), lines)]), header=FALSE, skip=1, col.names=header)
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
summary.bcftools.stats <- function(obj) {
    obj$SN
}


##' gplot.bcftools.stats
##'
##' Plot summary of bcftools stats results
##' @title plot.bcftools.stats
##' @param obj bcftools.stats object
##' @param which which plots to produce (defalt all)
##' @param ncol number of columns in grid plot
##' @param ... parameters passed to generic plot function
##' @return ggplot2 object, or a plot
##' @author Per Unneberg
gplot.bcftools.stats <- function(obj, which=c("SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP"), ncol=2, ...) {
    which <- match.arg(which, c("SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP"), several.ok=TRUE)
    message("Producing ", length(which), " plots")
    plist <- list()
    if ("SN" %in% which) {
        data <- obj$SN
        plist$SN <- ggplot(data=data, aes(y=key, x=value)) + geom_point() + ggtitle("Summary numbers") + theme_bw()
    }
    if ("TSTV" %in% which) {
        data <- gather(obj$TSTV, key, value, -id)
        tstv <- c("ts.tv", "ts.tv..1st.ALT.")
        plist$TSTV <- grid.arrange(
            ggplot(data=subset(data, key %in% tstv), aes(y=key, x=value)) + geom_point() + ggtitle("Ti/Tv") + theme_bw(),
            ggplot(data=subset(data, !(key %in% tstv)), aes(y=key, x=value)) + geom_point() + ggtitle("Ti/Tv, counts") + theme_bw(),
            nrow=2)
    }
    if ("SiS" %in% which) {
        data <- gather(obj$SiS, key, value, -id)
        plist$SiS <- ggplot(data=data, aes(y=key, x=value)) + geom_point() + ggtitle("Singleton stats") + theme_bw()
    }
    if ("AF" %in% which) {
        data <- gather(obj$AF, key, value, -allele.frequency, -id)
        plist$AF <- ggplot(data=data, aes(x=allele.frequency, y=value, color=key)) + geom_point() + ggtitle("Stats non-reference allele frequency") + theme_bw()
    }
    if ("QUAL" %in% which) {
        data <- gather(obj$QUAL, key, value, -Quality, -id)
        plist$QUAL <- ggplot(data=data, aes(x=Quality, y=value, color=key)) + geom_point() + ggtitle("Stats by quality") + theme_bw()
    }
    if ("IDD" %in% which) {
        data <- obj$IDD
        plist$IDD <- ggplot(data=data, aes(x=length..deletions.negative., y=count)) + geom_point() + ggtitle("InDel distribution") + xlab("Length (deletions negative)") + theme_bw()
    }
    if ("ST" %in% which) {
        data <- obj$ST
        plist$ST <- ggplot(data=data, aes(x=type, y=count)) + geom_point() + ggtitle("Substitution types") + theme_bw()
    }
    if ("DP" %in% which) {
        data <- gather(obj$DP, key, value, -bin, -id)
        plist$DP <- ggplot(data=data, aes(x=bin, y=value, color=key)) + geom_point() + ggtitle("Depth distribution") + theme_bw()
    }

    if (length(plist) > 1) {
        p <- grid.arrange(grobs=plist, ncol=ncol)
    } else {
        p <- plist[[1]]
    }
    p
}
