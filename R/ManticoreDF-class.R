.valid.ManticoreDF <- function(object) {
}

##' ManticoreDF
##'
##' Derived data frame class that requires that certain metadata
##' properties be present
##'
##'
##' @param object
##' @return
##' @author Per Unneberg
##'
##' @importFrom S4Vectors DataFrame
##'
setClass("ManticoreDF",
         contains = c("DataFrame"),
         validity = .valid.ManticoreDF,
         representation = representation(
             measurement.name = "character",
             application = "character"
         )
         )


setMethod("initialize", "ManticoreDF", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})


ManticoreDF <- function(..., measurement.name = character(), application = character()) {
    df <- DataFrame(...)
    x <- as(df, "ManticoreDF")
    x@measurement.name = measurement.name
    x@application = application
    x
}
