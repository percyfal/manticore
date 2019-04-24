fst <- system.file("extdata", "popoolation2", "p1_p2_w500.fst", package = "manticore")
fet <- system.file("extdata", "popoolation2", "p1_p2.fet", package = "manticore")
pwc <- system.file("extdata", "popoolation2", "p1_p2.sync_pwc", package = "manticore")
rc <- system.file("extdata", "popoolation2", "p1_p2.sync_rc", package = "manticore")
cmh <- system.file("extdata", "popoolation2", "p1_p2_p1_p2.cmh", package = "manticore")

context("test popoolation2 functions")

test_that("readPopoolation2 fet works", {
    message("reading readPopoolation2 fet")
    data <- readPopoolation2(fet)
})
test_that("readPopoolation2 fst works", {
    message("reading readPopoolation2 fst")
    data <- readPopoolation2(fst, "fst")
})
test_that("readPopoolation2 pwc works", {
    message("reading readPopoolation2 pwc")
    data <- readPopoolation2(pwc, "pwc")
})
test_that("readPopoolation2 rc works", {
    message("reading readPopoolation2 rc")
    data <- readPopoolation2(rc, "rc")
})
test_that("readPopoolation2 cmh works", {
    message("reading readPopoolation2 cmh")
    data <- readPopoolation2(cmh, "cmh")
})
