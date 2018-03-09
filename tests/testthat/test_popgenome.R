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
##     results.gr <- lapply(c("summary", stats), function(x) {getGenomeStats(scaffolds, x, quietly = TRUE, use.population.names = TRUE, out.format="GRanges")})
##     names(results.gr) <- c("summary", stats)

##     results.slide <- lapply(stats, function(x) {getGenomeStats(slides, x, quietly = TRUE, use.population.names = TRUE)})
##     names(results.slide) <- stats
##     results.slide.gr <- lapply(stats, function(x) {getGenomeStats(slides, x, quietly = TRUE, use.population.names = TRUE, out.format="GRanges")})
##     names(results.slide.gr) <- stats
##                                         ##save(slides, scaffolds, results, results.slide, scaffold2, scaffold13, slide2, slide13, populations.list, file="inst/extdata/popgenome.rda")
##     save(slides, scaffolds, scaffold2, scaffold13, slide2, slide13, populations.list, file="inst/extdata/popgenome.rda")
## }

if (all(unlist(lapply(c("PopGenome", "tidyr", "GenomicRanges", "tidyselect"),
                      requireNamespace, quietly = TRUE)))) {
    library("PopGenome")
    library("GenomicRanges")
    library("tidyr")
    library("tidyselect")
    fn <- system.file("extdata", "popgenome.rda", package = "nonmodelr")
    load(fn)
}

expect_ggplot <- function(data, gr) {
    p <- plot(data)
    eval(bquote(expect_is(p, "ggplot")))
    p <- boxplot(data)
    eval(bquote(expect_is(p, "ggplot")))
    p <- plot(data)
    eval(bquote(expect_is(p, "ggplot")))
    p <- boxplot(gr)
    eval(bquote(expect_is(p, "ggplot")))
}

expect_granges <- function(gr, col.names=sort(c("key", "value", "population")), r.width=c(340000, 10000)) {
    eval(bquote(expect_is(gr, "GRanges")))
    eval(bquote(expect_equal(width(gr)[1:2], r.width)))
    eval(bquote(expect_equal(sort(colnames(values(gr))), col.names)))
}

expect_pg <- function(data, type, col.names=sort(c("population", "ranges", "seqnames", "key", "value"))) {
    eval(bquote(expect_is(data, type)))
    eval(bquote(expect_equal(sort(colnames(data)), col.names)))
}

context("Test getGenomeStats functions")

test_that("summary data", {
    check_popgenome()
    col.names <- sort(c("n.sites", "n.biallelic.sites", "n.gaps",
                                              "n.unknowns", "n.valid.sites", "n.polyallelic.sites",
                        "trans.transv.ratio", "seqnames", "ranges", "population"))
    data <- getGenomeStats(scaffolds, "summary", out.format = "wide")
    expect_pg(data, "pg.summary", col.names = col.names)
    expect_equal(data$n.sites, c(340000, 10000))

    gr <- getGenomeStats(scaffolds, "summary", out.format = "GRanges")
    expect_granges(gr)
})

test_that("detail data", {
    check_popgenome()
    r.width <- rep(340000, 2)
    data <- getGenomeStats(scaffolds, "detail")
    expect_pg(data, "pg.detail")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("MDSD" %in% levels(as.factor(data$key)))

    gr <- getGenomeStats(scaffolds, "detail", out.format = "GRanges")
    expect_granges(gr, r.width = r.width)

    expect_ggplot(data, gr)
})

test_that("neutrality data", {
    check_popgenome()
    r.width <- rep(340000, 2)
    data <- getGenomeStats(scaffolds, "neutrality")
    expect_pg(data, "pg.neutrality")
    expect_true("pop 1" %in% levels(as.factor(data$population)))
    expect_true("Tajima.D" %in% levels(as.factor(data$key)))

    gr <- getGenomeStats(scaffolds, "neutrality", out.format = "GRanges")
    expect_granges(gr, r.width = r.width)

    expect_ggplot(data, gr)
})

test_that("fixed / shared", {
    check_popgenome()
    data.fixed <- getGenomeStats(scaffolds, "fixed")
    expect_pg(data.fixed, "pg.fixed")
    expect_true("pop1/pop2" %in% levels(as.factor(data.fixed$population)))

    gr.fixed <- getGenomeStats(scaffolds, "fixed", out.format = "GRanges")
    expect_granges(gr.fixed)
    expect_ggplot(data.fixed, gr.fixed)

    data.shared <- getGenomeStats(scaffolds, "shared")
    expect_pg(data.shared, "pg.shared")
    expect_true("pop1/pop2" %in% levels(as.factor(data.shared$population)))

    gr.shared <- getGenomeStats(scaffolds, "shared", out.format = "GRanges")
    expect_granges(gr.shared)

    expect_ggplot(data.shared, gr.shared)

    expect_false(identical(data.fixed, data.shared))
    expect_false(identical(gr.fixed, gr.shared))
})

test_that("sites", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "segregating.sites")
    expect_pg(data, "pg.segregating.sites")
    expect_true("pop 1" %in% levels(as.factor(data$population)))

    gr <- getGenomeStats(scaffolds, "segregating.sites", out.format = "GRanges")
    expect_granges(gr)

    expect_ggplot(data, gr)
})

test_that("diversity data", {
    check_popgenome()
    r.width <- rep(340000, 2)
    data <- getGenomeStats(scaffolds, "diversity")
    expect_pg(data, "pg.diversity")
    expect_true("nuc.diversity.within" %in% levels(as.factor(data$key)))
    expect_true("pop 1" %in% levels(as.factor(data$population)))

    gr <- getGenomeStats(scaffolds, "diversity", out.format = "GRanges")
    expect_granges(gr, r.width = r.width)

    expect_ggplot(data, gr)
})

test_that("diversity.between data", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "diversity.between")
    expect_pg(data, "pg.diversity.between")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))

    gr <- getGenomeStats(scaffolds, "diversity.between", out.format = "GRanges")
    expect_granges(gr)

    expect_ggplot(data, gr)
})

test_that("F_ST", {
    check_popgenome()
    data <- getGenomeStats(scaffolds, "F_ST")
    expect_pg(data, "pg.F_ST")

    gr <- getGenomeStats(scaffolds, "F_ST", out.format = "GRanges")
    expect_granges(gr)

    expect_ggplot(data, gr)
})

test_that("F_ST.pairwise", {
    check_popgenome()
    r.width <- rep(340000, 2)
    data <- getGenomeStats(scaffolds, "F_ST.pairwise")
    expect_pg(data, "pg.F_ST.pairwise")
    expect_true("pop1/pop2" %in% levels(as.factor(data$population)))

    gr <- getGenomeStats(scaffolds, "F_ST.pairwise", out.format = "GRanges")
    expect_granges(gr, r.width = r.width)

    expect_ggplot(data, gr)
})
