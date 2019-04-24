context("Windows tests")

x <- Windows(
    seqnames = c(rep("1", 2), rep("2", 3)),
    IRanges::IRanges(start = c(seq(1, 20, 10), seq(11, 40, 10)),
                     end = c(seq(10, 20, 10), seq(20, 40, 10))))

test_that("making coordinates preserves shape", {
    wc <- makeCoordinates(x)
    expect_equal(length(wc), length(x))
})
