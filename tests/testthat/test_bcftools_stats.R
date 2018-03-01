context("Test parsing bcftools output, S4")
fn <- file.path("../../inst/extdata/medium.call.stats")
S4stats <- readBcftoolsStats(fn)

test_that("new S4 bcftools instance", {
    expect_equal(class(S4stats)[1], "BcftoolsStats")
    expect_equal(S4stats@label, fn)
})


context("Test parsing bcftools output, S3")
S3stats <- read.bcftools.stats(fn)
test_that("new S3 bcftools instance", {
    expect_equal(class(S3stats)[1], "bcftools.stats")
    expect_equal(S3stats$label, fn)
})
