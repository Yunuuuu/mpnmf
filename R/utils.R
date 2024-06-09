`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

is_scalar <- function(x) length(x) == 1L
