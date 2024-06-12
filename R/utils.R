`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

is_scalar <- function(x) length(x) == 1L

# deal with RcppML `nmf` object, since `nmf` is not exported in CRAN version
is_RcppML_nmf <- function(x) {
    methods::is(x, "nmf") && identical(attr(class(x), "package"), "RcppML")
}

basis <- function(x) {
    if (is_RcppML_nmf(x)) {
        x@w
    } else {
        NMF::basis(x)
    }
}

nbasis <- function(x) {
    if (is_RcppML_nmf(x)) {
        ncol(x@w)
    } else {
        NMF::nbasis(x)
    }
}

hasBasis <- function(x) nrow(basis(x)) && nbasis(x)

coef <- function(x) {
    if (is_RcppML_nmf(x)) {
        x@h
    } else {
        NMF::coef(x)
    }
}

hasCoef <- function(x) nbasis(x) && ncol(coef(x))
