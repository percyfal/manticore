overlapByWindows <- function(windows, subject, ignore.strand = TRUE, ...) {
    message(paste0("looking for overlaps in ", length(windows), " windows"))
    hits <- findOverlaps(windows, subject, ignore.strand = ignore.strand)
    message(paste0("found ", length(hits), " overlaps"))
    x.start <- apply(
        cbind(
            start(subject[subjectHits(hits)]),
            start(windows[queryHits(hits)])),
        1, max)
    x.end <- apply(
        cbind(
            end(subject[subjectHits(hits)]),
            end(windows[queryHits(hits)])),
        1, min)

    x <- reduce(GRanges(queryHits(hits), IRanges(start = x.start, end = x.end)))
    olap <- tapply(width(x), seqnames(x), sum)
    y <- rep(0, length(windows))
    y[as.integer(names(olap))] <- olap
    y
}
