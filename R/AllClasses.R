##' Collect Genome-wide statistics
##'
##' Subclass of RangedSummarizedExperiment
##'
##' @export
##' @rdname GStats-class
##'
setClass("GStats",
         contains = c("RangedSummarizedExperiment"),
         representation = representation(
             statistics = "character"
         )
         )

.valid.GStats <- function(x) {
    if (length(x)) {
        if (!identical(sort(colnames(colData(x))), sort(c("population", "statistic"))))
            return("colData must be indexed by factors 'population' and 'statistic'")
    }
    NULL
}
setValidity2("GStats", .valid.GStats)

##' List of GENOME instances
##'
##' Subclass of S4Vectors SimpleList, where each entry is an PopGenome
##' GENOME instance
##'
##' @export
##' @rdname GENOMEList-class
##'
setClass("GENOMEList",
         contains = "SimpleList",
         prototype = prototype(elementType = "GENOME")
         )
