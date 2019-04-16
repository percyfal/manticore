dmel.df <- do.call("rbind", apply(expand.grid(c("A", "B"), c("D", "pi", "theta")), 1, function(x){data.frame(x[1], x[2], system.file("extdata", paste0("popoolation/dmel.", x[1], ".", x[2], ".txt.gz"), package="manticore"))}))
colnames(dmel.df) <- c("sample", "measure", "filename")

context("test popoolation functions")

test_that("readVarianceSliding works", {
    dmel <- as.list(as.character(subset(dmel.df, sample == "A")$filename))
    names(dmel) <- subset(dmel.df, sample == "A")$measure
    tmp <- lapply(names(dmel), function(x) {readVarianceSliding(dmel[[x]], x, sample = "A")})
})


test_that("VarianceSlidingAssays works", {
    tmp <- VarianceSlidingAssays(dmel.df)
})
