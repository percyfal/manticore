##' makeCoordinates
##'
##'
##'
##' @param x Windows object
##' @param mapping data frame that orders seqnames and possibly maps
##'     them to a chromosome name
##' @return updated Windows object with index and pos columns
##' @author Per Unneberg
##'
##' @importFrom IRanges slidingWindows
##' @importFrom GenomeInfoDb seqinfo seqnames seqlengths
##'
makeCoordinates <- function(x, mapping = NULL) {
    message("Making coordinates")
    w0 <- width(x)[1]
    if (any(is.na(seqlengths(x)))) {
        slen <- tapply(end(x), seqnames(x), max)
    } else {
        slen <- seqlengths(x)
    }
    y <- GenomicRanges::GRanges(seqnames = levels(seqnames(x)),
                                ranges = IRanges::IRanges(start = 1, end = slen))
    w <- unlist(IRanges::slidingWindows(y, width = w0, step = w0))
    mcols(x)$index <- match(x, w)
    mcols(x)$pos <- cumsum(width(w))[mcols(x)$index] -w0 + 1
    x
}

##' @title addSeqnamesColor
##'
##' @description Add color mapping to seqnames
##'
##' @param x
##' @return updated Windows object with color column
##' @author Per Unneberg
##'
addSeqnamesColor <- function(x, mapping = NULL, n.levels = 2) {
    if ("colour" %in% names(x))
        return(x)
    if (is.null(mapping)) {
        mcols(x)$colour <- (as.integer(seqnames(x)) + n.levels - 1) %% n.levels + 1
    } else {
        mcols(x)$colour <- (as.integer(seqnames(x)) + n.levels - 1) %% n.levels + 1
    }
    mcols(x)$colour <- as.factor(mcols(x)$colour)
    return(x)
}
