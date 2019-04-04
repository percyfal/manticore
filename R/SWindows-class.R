### ======================================================================
### SWindows class
### Segregating windows with information on segregating and good sites
### ----------------------------------------------------------------------
###


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
             coverage = "numeric",
             sites = "integer",
             segregating.sites = "integer"
         ),
         prototype = prototype(
             coverage = numeric(),
             sites = integer(),
             segregating.sites = integer()
         )
         )

setMethod("parallelSlotNames", "SWindows",
          function(x) c("rowRanges", callNextMethod())
          )

## setMethod(GenomicRanges:::extraColumnSlotNames, "SWindows",
##           function(x) {
##     c("coverage", "sites", "segregating.sites")
## })

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
    new("SWindows", w, coverage = coverage, sites = sites, segregating.sites = segregating.sites)
}


### ----------------------------------------------------------------------
### Functions
###

#sites <- function
