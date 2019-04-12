### ======================================================================
### SWindows class
### Segregating windows with information on segregating and good sites
### ----------------------------------------------------------------------
###

setClassUnion("integer_OR_NA", c("integer", "logical"))
setClassUnion("numeric_OR_NA", c("numeric", "logical"))

##' SWindows class
##'
##' Class for windows carrying additional parallel slots with
##' information about number of good sites, window coverage and
##' segregating sites
##'
##' @export
##' @rdname SWindows-class
##'
##' @author Per Unneberg
##'
setClass("SWindows",
         contains = c("Windows"),
         representation = representation(
             coverage = "numeric_OR_NA",
             sites = "integer_OR_NA",
             segregating.sites = "integer_OR_NA"
         ),
         prototype = prototype(
             coverage = numeric(),
             sites = integer(),
             segregating.sites = integer()
         )
         )


##' parallelSlotNames
##'
##'
##'
##' @param x SWindows object
##'
##' @author Per Unneberg
##' @importFrom S4Vectors parallelSlotNames
##'
setMethod("parallelSlotNames", "SWindows",
          function(x) c("coverage", "sites", "segregating.sites", callNextMethod())
          )


setMethod(GenomicRanges:::extraColumnSlotNames, "SWindows",
          function(x) {
    c("coverage", "sites", "segregating.sites")
})

### ----------------------------------------------------------------------
### Constructors
###


##' SWindows constructor
##'
##' Create an SWindows object
##'
##' @param ... arguments to pass on to parent class Windows
##' @param coverage coverage per window
##' @param sites number of good sites per window
##' @param segregating.sites number of segregating sites per window
##'
##' @rdname SWindows-class
##' @return SWindows class
##' @author Per Unneberg
##'
SWindows <- function(..., coverage = numeric(),
                     sites = integer(),
                     segregating.sites = integer())
{
    w <- Windows(...)
    if (missing(coverage))
        coverage <- rep(NA, length(w))
    if (missing(sites))
        sites <- rep(NA, length(w))
    if (missing(segregating.sites))
        segregating.sites <- rep(NA, length(w))
    new("SWindows", w, coverage = coverage, sites = sites, segregating.sites = segregating.sites)
}


### ----------------------------------------------------------------------
### Functions
###
