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

if (all(unlist(lapply(c("PopGenome", "tidyr", "GenomicRanges"), requireNamespace, quietly = TRUE)))) {
    library("PopGenome")
    library("GenomicRanges")
    library("tidyr")
    ## TODO: prepare this data in a data file
    populations.list <- list(
        CHS = c("CHS.HG00512", "CHS.HG00513"),
        PUR = c("PUR.HG00731", "PUR.HG00733"),
        YRI = c("YRI.NA19238", "YRI.NA19239")
    )
    fn <- system.file("extdata", "medium.call.biallelic.vcf.gz", package="nonmodelr")
    scaffold2 <- readVCF(fn, frompos = 1, topos = 340000,
                         numcols = 1000, tid = "scaffold2")
    scaffold13 <- readVCF(fn, frompos = 1, topos = 10000,
                          numcols = 1000, tid = "scaffold13")
    scaffold2@region.names <- "scaffold2"
    scaffold13@region.names <- "scaffold13"
    scaffold2 <- set.populations(scaffold2, populations.list, diploid = TRUE)
    scaffold13 <- set.populations(scaffold13, populations.list, diploid = TRUE)
    scaffold2 <- genomewide.stats(scaffold2)
    scaffold13 <- genomewide.stats(scaffold13)
    scaffolds <- list(scaffold2 = scaffold2, scaffold13 = scaffold13)

    slide2 <- sliding.window.transform(scaffold2,
                                       10000,
                                       10000,
                                       type = 2, whole.data = FALSE)
    slide13 <- sliding.window.transform(scaffold13,
                                        1000,
                                        1000,
                                        type = 2, whole.data = FALSE)
    slide2 <- genomewide.stats(slide2)
    slide13 <- genomewide.stats(slide13)
    slides <- list(slide2 = slide2, slide13 = slide13)
}

context("Test getGenomeStats functions")

test_that("summary data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "summary", out.format = "wide")
    expect_is(data, "pg.summary")
    expect_equal(data$n.sites, c(340000, 10000))
    expect_equal(sort(colnames(data)), sort(c("n.sites", "n.biallelic.sites", "n.gaps",
                                              "n.unknowns", "n.valid.sites", "n.polyallelic.sites",
                                              "trans.transv.ratio", "seqnames")))
})

test_that("detail data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "detail")
    expect_is(data, "pg.detail")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("MDSD" %in% levels(as.factor(data$key)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
})

test_that("neutrality data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "neutrality")
    expect_is(data, "pg.neutrality")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("Tajima.D" %in% levels(as.factor(data$key)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
})

test_that("fixed / shared", {
    check_popgenome()
    data.fixed <- getGenomeStats(scaffolds, "fixed")
    expect_is(data.fixed, "pg.fixed")
    expect_true("pop1/pop2" %in% levels(as.factor(data.fixed$population)))
    expect_equal(sort(colnames(data.fixed)),
                 sort(c("population", "ranges", "seqnames", "value")))

    data.shared <- getGenomeStats(scaffolds, "shared")
    expect_is(data.shared, "pg.shared")
    expect_true("pop1/pop2" %in% levels(as.factor(data.shared$population)))
    expect_equal(sort(colnames(data.shared)),
                 sort(c("population", "ranges", "seqnames", "value")))

    expect_false(identical(data.fixed, data.shared))

})

test_that("sites", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "segregating.sites")
    expect_is(data, "pg.segregating.sites")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value")))
})

test_that("diversity data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "diversity")
    expect_is(data, "pg.diversity")
    expect_true("nuc.diversity.within" %in% levels(as.factor(data$key)))
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "key", "value")))
})

test_that("diversity.between data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "diversity.between")
    expect_is(data, "pg.diversity.between")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value")))
})

test_that("F_ST", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "F_ST")
    expect_is(data, "pg.F_ST")
    expect_true(!("population" %in% colnames(data)))
    expect_equal(sort(colnames(data)),
                 sort(c("ranges", "seqnames", "value", "key")))
})

test_that("F_ST.pairwise", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "F_ST.pairwise")
    expect_is(data, "pg.F_ST.pairwise")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))
    expect_equal(sort(colnames(data)),
                 sort(c("population", "ranges", "seqnames", "value", "key")))
})
