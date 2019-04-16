a <- unlist(slidingWindows(GRanges(seqnames="foo", ranges=IRanges(start=1, end=10000)), 1000, 1000))
j.start <- c(428, 503, 1102, 1115, 1579, 2993, 4815, 7466, 7929, 7968)
j.end <- c(681, 588, 1305, 1222, 1634, 3678, 4952, 8158, 8687, 8298)
b <- sort(GRanges(seqnames="foo", ranges=IRanges(start=j.start, end=j.end)))

test_that("windows overlap works", {
    olap <- overlapByWindows(a, b)
    expect_equal(olap[1], width(reduce(b))[1])
    expect_equal(olap[3], 8)
    expect_equal(olap[6], 0)
})
