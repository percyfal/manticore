##' Collect Genome-wide statistics
##'
##' Subclass of RangedSummarizedExperiment
##'
##' @export
##' @rdname GStats-class
##'
##' @import SummarizedExperiment
##' @importFrom methods validObject
##'
setClass("GStats",
         contains = c("RangedSummarizedExperiment"),
         representation = representation(
             statistics = "character",
             application = "character"
         )
         )

.valid.GStats <- function(object) {
    if (length(object)) {
        if (!identical(sort(colnames(colData(object))), sort(c("population", "statistic"))))
            return("colData must be indexed by factors 'population' and 'statistic'")
        if (!identical(assayNames(object), c("data")))
            return("only one assay ('data') allowed")
    }
    NULL
}
setValidity("GStats", .valid.GStats)

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
