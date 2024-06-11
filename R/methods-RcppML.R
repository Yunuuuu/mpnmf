#' Methods for RcppML nmf
#'
#' Add [NMF][NMF] package methods for RcppML [nmf][RcppML::nmf] object.
#'
#' @param object A [nmf][RcppML::nmf] object.
#' @param ... Not used currently.
#' @seealso
#' - [basis][NMF::basis]
#' - [coef][NMF::coef]
#' @name method-nmf
NULL

#' @importClassesFrom RcppML nmf
#' @importFrom NMF basis
#' @export
#' @rdname method-nmf
methods::setMethod("basis", "nmf", function(object, ...) object@w)

#' @importFrom NMF coef
#' @export
#' @rdname method-nmf
methods::setMethod("coef", "nmf", function(object, ...) object@h)
