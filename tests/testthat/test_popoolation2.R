fst <- system.file("extdata", "popoolation2", "p1_p2_w500.fst", package = "manticore")
fet <- system.file("extdata", "popoolation2", "p1_p2.fet", package = "manticore")
pwc <- system.file("extdata", "popoolation2", "p1_p2.sync_pwc", package = "manticore")
rc <- system.file("extdata", "popoolation2", "p1_p2.sync_rc", package = "manticore")
cmh <- system.file("extdata", "popoolation2", "p1_p2_p1_p2.cmh", package = "manticore")

context("test popoolation2 functions")

test_that("readPopoolation2 works", {
    message("reading readPopoolation2")
    data <- readPopoolation2(fst, assay)
    head(data)
})
