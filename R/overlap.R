## a <- unlist(slidingWindows(GRanges(seqnames="foo", ranges=IRanges(start=1, end=100000)), 1000, 1000))
## j.start <- sample(max(end(a)), 100)
## j.end <- j.start + sample(1000, 100)
## b <- sort(GRanges(seqnames="foo", ranges=IRanges(start=j.start, end=j.end)))
overlapByWindows <- function(windows, subject, ignore.strand = TRUE, ...) {
    hits <- findOverlaps(windows, subject, ignore.strand = ignore.strand)
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
