### ======================================================================
### Windows class
### ----------------------------------------------------------------------
###


##' Windows representation
##'
##' Windows object to represent windowed analyses. Contains a slot for
##' window size. Note that edge effects are allowed such that the last
##' windows may be smaller than the actual window size.
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
    if (!length(window.size)) {
        window.size <- start(gr)[2] - start(gr)[1]
    }
    start(gr) <- start(gr) - window.size / 2 + 1
    end(gr) <- end(gr) + window.size / 2
    gr <- trim(gr)
    new("Windows", gr, window.size = as.integer(window.size))
}
