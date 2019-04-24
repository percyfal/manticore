check_popgenome <- function() {
    if (!requireNamespace("PopGenome", quietly = TRUE)) {
        skip("PopGenome package not available")
    }
}

if (requireNamespace("PopGenome", quietly = TRUE))
    library("PopGenome")

fn <- system.file("extdata", "popgenome.rda", package = "manticore")
load(fn)

context("Test PopGenome methods")

test_that("parsing slide data is ok", {
    check_popgenome()
    wselist <- PopGenome(slides)

})
