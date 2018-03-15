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
