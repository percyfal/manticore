dmel.pi <- system.file("extdata", "popoolation/dmel.pi.txt.gz", package = "manticore")
dmel.D <- system.file("extdata", "popoolation/dmel.D.txt.gz", package = "manticore")
dmel.theta <- system.file("extdata", "popoolation/dmel.theta.txt.gz", package = "manticore")
dmel <- list(dmel.pi, dmel.D, dmel.theta)
names(dmel) <- c("pi", "D", "theta")
dmel.df <- data.frame(filename = rep(as.character(dmel), 2),
                      sample = rep(c("A", "B"), each=3),
                      measure = rep(as.character(names(dmel)), 2))

context("test popoolation functions")

test_that("readVarianceSliding works", {
    tmp <- lapply(names(dmel), function(x) {readVarianceSliding(dmel[[x]], x)})
})


test_that("VarianceSlidingAssays works", {
    tmp <- VarianceSlidingAssays(dmel.df)
})
