check_popgenome <- function() {
    if (!requireNamespace("PopGenome", quietly = TRUE)) {
        skip("PopGenome package not available")
    }
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        skip("GenomicRanges package not available")
    }
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        skip("tidyr package not available")
    }

}

## save.popgenome <- FALSE
## if (save.popgenome) {
##     ## Code to save popgenome.rda
##     populations.list <- list(
##         CHS = c("CHS.HG00512", "CHS.HG00513"),
##         PUR = c("PUR.HG00731", "PUR.HG00733"),
##         YRI = c("YRI.NA19238", "YRI.NA19239")
##     )
##     fn <- system.file("extdata", "medium.call.biallelic.vcf.gz", package="nonmodelr")
##     scaffold2 <- readVCF(fn, frompos = 1, topos = 340000,
##                          numcols = 1000, tid = "scaffold2")
##     scaffold13 <- readVCF(fn, frompos = 1, topos = 10000,
##                           numcols = 1000, tid = "scaffold13")
##     scaffold2@region.names <- "scaffold2"
##     scaffold13@region.names <- "scaffold13"
##     scaffold2 <- set.populations(scaffold2, populations.list, diploid = TRUE)
##     scaffold13 <- set.populations(scaffold13, populations.list, diploid = TRUE)
##     scaffold2 <- genomewide.stats(scaffold2)
##     scaffold13 <- genomewide.stats(scaffold13)
##     scaffolds <- list(scaffold2 = scaffold2, scaffold13 = scaffold13)

##     slide2 <- sliding.window.transform(scaffold2,
##                                        10000,
##                                        10000,
##                                        type = 2, whole.data = FALSE)
##     slide13 <- sliding.window.transform(scaffold13,
##                                         1000,
##                                         1000,
##                                         type = 2, whole.data = FALSE)
##     slide2 <- genomewide.stats(slide2)
##     slide13 <- genomewide.stats(slide13)
##     slides <- list(scaffold2 = slide2, scaffold13 = slide13)
##     stats <- c("detail", "neutrality", "fixed", "shared", "diversity", "diversity.between",
##                "F_ST", "F_ST.pairwise", "segregating.sites")
##     results <- lapply(c("summary", stats), function(x) {getGenomeStats(scaffolds, x, quietly = TRUE, use.population.names = TRUE)})
##     names(results) <- c("summary", stats)
##     results.slide <- lapply(stats, function(x) {getGenomeStats(slides, x, quietly = TRUE, use.population.names = TRUE)})
##     names(results.slide) <- stats
##     save(slides, scaffolds, results, results.slide, scaffold2, scaffold13, slide2, slide13, populations.list, file="inst/extdata/popgenome.rda")
## }

if (all(unlist(lapply(c("PopGenome", "tidyr", "GenomicRanges"),
                      requireNamespace, quietly = TRUE)))) {
    library("PopGenome")
    library("GenomicRanges")
    library("tidyr")
    fn <- system.file("extdata", "popgenome.rda", package = "nonmodelr")
    load(fn)
}

context("Test getGenomeStats functions")

test_that("summary data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "summary", out.format = "wide")
    expect_is(data, "pg.summary")
    expect_equal(data$n.sites, c(340000, 10000))
    expect_equal(sort(colnames(data)), sort(c("n.sites", "n.biallelic.sites", "n.gaps",
                                              "n.unknowns", "n.valid.sites", "n.polyallelic.sites",
                                              "trans.transv.ratio", "seqnames", "ranges", "population")))
    gr <- getGenomeStats(scaffolds, "summary", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("detail data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "detail")
    expect_is(data, "pg.detail")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("MDSD" %in% levels(as.factor(data$key)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
    gr <- getGenomeStats(scaffolds, "detail", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], rep(340000, 2))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("neutrality data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "neutrality")
    expect_is(data, "pg.neutrality")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("Tajima.D" %in% levels(as.factor(data$key)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
    gr <- getGenomeStats(scaffolds, "neutrality", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], rep(340000, 2))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("fixed / shared", {
    check_popgenome()
    data.fixed <- getGenomeStats(scaffolds, "fixed")
    expect_is(data.fixed, "pg.fixed")
    expect_true("pop1/pop2" %in% levels(as.factor(data.fixed$population)))
    expect_equal(sort(colnames(data.fixed)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))
    gr <- getGenomeStats(scaffolds, "fixed", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))


    data.shared <- getGenomeStats(scaffolds, "shared")
    expect_is(data.shared, "pg.shared")
    expect_true("pop1/pop2" %in% levels(as.factor(data.shared$population)))
    expect_equal(sort(colnames(data.shared)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))

    expect_false(identical(data.fixed, data.shared))
    gr.shared <- getGenomeStats(scaffolds, "shared", out.format = "GRanges")
    expect_is(gr.shared, "GRanges")
    expect_equal(width(gr.shared)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr.shared))),
                 sort(c("key", "value", "population")))


})

test_that("sites", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "segregating.sites")
    expect_is(data, "pg.segregating.sites")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))
    gr <- getGenomeStats(scaffolds, "segregating.sites", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))

})

test_that("diversity data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "diversity")
    expect_is(data, "pg.diversity")
    expect_true("nuc.diversity.within" %in% levels(as.factor(data$key)))
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
    gr <- getGenomeStats(scaffolds, "diversity", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], rep(340000, 2))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("diversity.between data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "diversity.between")
    expect_is(data, "pg.diversity.between")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))
    gr <- getGenomeStats(scaffolds, "diversity.between", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("F_ST", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "F_ST")
    expect_is(data, "pg.F_ST")
    expect_equal(sort(colnames(data)),
                 sort(c("ranges", "seqnames", "value", "key", "population")))
    gr <- getGenomeStats(scaffolds, "F_ST", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], c(340000, 10000))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})

test_that("F_ST.pairwise", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "F_ST.pairwise")
    expect_is(data, "pg.F_ST.pairwise")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))
    gr <- getGenomeStats(scaffolds, "F_ST.pairwise", out.format = "GRanges")
    expect_is(gr, "GRanges")
    expect_equal(width(gr)[1:2], rep(340000, 2))
    expect_equal(sort(colnames(values(gr))),
                 sort(c("key", "value", "population")))
})
