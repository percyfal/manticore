if (!requireNamespace("PopGenome", quietly = TRUE)) {
    skip("PopGenome package not available")
}


if (requireNamespace("PopGenome", quietly = TRUE))
    library("PopGenome")

fn <- system.file("extdata", "popgenome.rda", package = "nonmodelr")
load(fn)
scaffolds.gl <- GENOMEList(scaffolds)

gs <- GStats(scaffolds.gl)

context("GStats functionality")

test_that("asGRanges works", {
    gr <- asGRanges(gs)
    expect_is(gr, "GRanges")
})
