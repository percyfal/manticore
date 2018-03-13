##' identify_outliers
##'
##' ientify outliers in a data frame
##' @title identify_outliers
##' @param data data frame
##' @param formula formula to fit
##' @param key factor over which to perform fits
##' @param method fitting method
##' @param level level of confidence interval to use
##' @param ... additional arguments
##' @return data frame with fitted model values, se, residuals and key
##' @author Per Unneberg
##' @export
identify_outliers <- function(data, formula, key=NULL, method="loess", level=0.99, ...) {
    ## See https://stackoverflow.com/questions/33082901/find-points-over-and-under-the-confidence-interval-when-using-geom-stat-geom-s for loess solution
    method <- match.arg(method, c("loess", "lm"))
    f <- match.fun(method)
    if (!is.null(key)) {
        conf <- do.call("rbind", lapply(levels(factor(data[[key]])),
                                        function(ll) {
                                     x <- data[data[[key]] == ll, ];
                                     fit <- f(formula, x);
                                     y <- predict(fit, se = TRUE);
                                     se.fit <- y$se.fit * qt(level / 2 + .5, y$df);
                                     data.frame(fit = y$fit, se.fit = se.fit, upper = y$fit + se.fit, lower = y$fit - se.fit, residuals = residuals(fit), key = ll)}))
    } else {
        fit <- f(value ~ seqlengths, data);
        y <- predict(fit, se = TRUE);
        se.fit <- y$se.fit * qt(level / 2 + .5, y$df);
        conf <- data.frame(fit = y$fit, se.fit = se.fit, upper = y$fit + se.fit, lower = y$fit - se.fit, residuals = residuals(fit), key = NA);
    }
    conf
}
