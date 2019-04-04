### ======================================================================
### Windows class
### ----------------------------------------------------------------------
###


##' Windows representation
##'
##' Windows object to represent windowed analyses. Contains a slot for
##' window size.
##'
##' @export
##' @rdname Windows-class
##'
##' @author Per Unneberg
##'
setClass("Windows",
         contains = c("GRanges"),
         representation = representation(
             window.size = "integer"
         ),
         prototype = prototype(
             window.size = integer()
         )
         )

### ----------------------------------------------------------------------
### Constructors
###

##' Windows
##'
##' Create Windows class
##'
##' @param ... parameters to be passed on to parent class GRanges
##' @param window.size window size of analysis
##'
##' @rdname Windows-class
##' @return Windows object
##' @author Per Unneberg
##'
Windows <- function(..., window.size = integer())
{
    gr <- GRanges(...)
    new("Windows", gr, window.size = window.size)
}
