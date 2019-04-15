w <- Windows(seqnames = c("1", "2"), IRanges(start = c(1, 10), end = c(2,20)))

test_that("making coordinates preserves shape", {
    wc <- makeCoordinates(w)
    message(wc)
})
