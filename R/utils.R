
## Defaults for NULL values
`%||%` <- function(a, b) if (is.null(a)) b else a

## Remove NULLs from a list
compact <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}
