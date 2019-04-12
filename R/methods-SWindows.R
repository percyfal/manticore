##' sites
##'
##' Get the number of sites
##'
##' @param obj SWindows object
##' @return Vector
##'
##' @author Per Unneberg
##' @export
##'
setMethod("sites", "SWindows", function(obj) {
    x <- integer()
    if (!all(is.na(obj@sites)))
        x <- obj@sites
    else if (!all(is.na(obj@coverage)))
        x <- as.integer(obj@coverage * width(obj))
    else
        NULL
    x
})
