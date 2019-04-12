## Return Windows with additional columns pos, index and chr
## Mapping is from seqnames to chr in case pseudomapping present
makeCoordinates <- function(x, mapping = NULL) {
    w0 <- width(x)[1]
    mcols(x)$pos <- cumsum(width(x)) - w0 + width(x) / 2
    x
}
