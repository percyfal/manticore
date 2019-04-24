##' @rdname overlapByWindows
##' @title overlapByWindows
##'
##' @description find overlaps by windows
##'
##' @param windows windows object
##' @param subject subject
##' @param ignore.strand whether or not to ignore strand
##'
##' @return
##' @author Per Unneberg
##'
##'
overlapByWindows <- function(windows, subject, ignore.strand = TRUE) {
    message(paste0("looking for overlaps in ", length(windows), " windows"))
    hits <- IRanges::findOverlaps(windows, subject, ignore.strand = ignore.strand)
    message(paste0("found ", length(hits), " overlaps"))
    x.start <- apply(
        cbind(
            start(subject[S4Vectors::subjectHits(hits)]),
            start(windows[S4Vectors::queryHits(hits)])),
        1, max)
    x.end <- apply(
        cbind(
            end(subject[S4Vectors::subjectHits(hits)]),
            end(windows[S4Vectors::queryHits(hits)])),
        1, min)

    x <- IRanges::reduce(GRanges(S4Vectors::queryHits(hits), IRanges::IRanges(start = x.start, end = x.end)))
    olap <- tapply(width(x), seqnames(x), sum)
    y <- rep(0, length(windows))
    y[as.integer(names(olap))] <- olap
    y
}
