check_popgenome <- function() {
    if (!requireNamespace("PopGenome", quietly = TRUE)) {
        skip("PopGenome package not available")
    }
}

if (requireNamespace("PopGenome", quietly = TRUE))
    library("PopGenome")

fn <- system.file("extdata", "popgenome.rda", package = "manticore")
load(fn)


## expect_ggplot <- function(gs) {
##     p <- gplot(gs)
##     eval(bquote(expect_is(p, "ggplot")))
##     p <- gboxplot(gs)
##     eval(bquote(expect_is(p, "ggplot")))
## }

## expect_gstats <- function(gs, r.width = c(340000, 10000)) {
##     eval(bquote(expect_is(gs, "GStats")))
##     eval(bquote(expect_equal(width(gs)[1:2], r.width)))
## }

context("Test PopGenome methods")

## test_that("getting summary data ok", {
##     check_popgenome()
##     col.names <- sort(c("n.sites", "n.biallelic.sites", "n.gaps",
##                         "n.unknowns", "n.valid.sites", "n.polyallelic.sites",
##                         "trans.transv.ratio", "seqnames", "ranges", "population"))

##     gs <- GStats(scaffolds, statistics = "summary")
##     expect_gstats(gs)
## })

## test_that("getting detail data ok", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "detail")
##     expect_gstats(gs)
## })

## test_that("plotting detail data ok", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "detail")
##     expect_ggplot(gs)
## })

## test_that("getting neutrality data ok", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "neutrality")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })

## test_that("getting fixed / shared data ok", {
##     check_popgenome()
##     gs.fixed <- GStats(scaffolds, statistics = "fixed")
##     expect_gstats(gs.fixed)
##     ## expect_ggplot(gs.fixed)

##     gs.shared <- GStats(scaffolds, statistics = "shared")
##     expect_gstats(gs.shared)
##     ## expect_ggplot(gs.shared)

##     expect_false(identical(gs.fixed, gs.shared))
## })

## test_that("getting segregating sites ok", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "segregating.sites")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })

## test_that("getting diversity data ok", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "diversity")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })

## test_that("getting diversity.between data works", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "diversity.between")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })

## test_that("getting F_ST data works", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "F_ST")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })

## test_that("getting F_ST.pairwise data works", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "F_ST.pairwise")
##     expect_gstats(gs)
##     ## expect_ggplot(gs)
## })


## test_that("coldata order corresponds to variable order", {
##     check_popgenome()
##     gs <- GStats(slides, statistics = "detail")
##     cdata <- paste(gs$population, gs$statistic, sep = "_")
##     expect_equal(rownames(colData(gs)), cdata)
## })

## test_that("GStats object with wrong slots fail", {
##     check_popgenome()
##     gs <- GStats(scaffolds, statistics = "detail")
##     expect_error(gs@foo <- "bar")
## })

## test_that("GStats with GRanges argument works", {
##     check_popgenome()
##     gr <- GRanges(seqnames=c("scaffold2"), ranges=IRanges(start=1e5, end=2e5))
##     gs <- GStats(scaffolds, gr)
##     expect_equal(rowRanges(gs)$sites[1], width(gr))
## })
