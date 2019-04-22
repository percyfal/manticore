##' sites
##'
##' @export
##'
##' @describeIn sites Get sites for an SWindows object
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

##' @rdname window.size
##' @description Get window.size from a Windows object
##'
##' @param obj Windows object
##'
setMethod("window.size", "SWindows",
          function(obj) {
    return (obj@window.size)
})
