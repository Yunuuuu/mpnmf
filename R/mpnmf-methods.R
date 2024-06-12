#' Methods for `mpnmf` object.
#' 
#' `mp_signatures`: A list of features to define the meta program signature.
#'
#' `mp_samples`: A list of samples indicates where the meta program is from.
#'
#' `mp_coverage`: A numeric indicates the fraction of samples the meta program
#'                was detected.
#'
#' `mp_identity`: A list or an atomic character indicates  which program the
#'                cells were allocated.
#' 
#' @param x A [mpnmf][mp] object.
#' @param s_min A scalar integer indicates the minimal number of samples or a
#' scalar numeric (`0 < s_min < 1`) indicates the minimal proportion of samples
#' to define the meta program. If `NULL`, no filters will be applied. Default:
#' `1/3`.
#' @name mpnmf-method
NULL

#' @param ... Not used currently.
#' @export
#' @rdname mpnmf-method
print.mpnmf <- function(x, s_min = 1 / 3, ...) {
    print(mp_programs(x, s_min))
}

#' @export
#' @rdname mpnmf-method
mp_programs <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_programs[index]
}

#' @param flatten A boolean value indicates whether to return an atomic
#' characters instead of a list.
#' @export
#' @rdname mpnmf-method
mp_identity <- function(x, s_min = 1 / 3, flatten = FALSE) {
    mp_programs <- mp_programs(x, s_min)
    assert_bool(flatten)
    identity_to_mp <- structure(
        rep(names(mp_programs), times = lengths(mp_programs)),
        names = unlist(mp_programs, FALSE, FALSE)
    )
    identity <- attr(x, "identity")
    out <- lapply(identity, function(index) {
        out <- identity_to_mp[index]
        if (!is.null(names(index))) names(out) <- names(index)
        out
    })
    if (flatten) {
        names(out) <- NULL
        out <- unlist(out, recursive = FALSE, use.names = TRUE)
    }
    out
}

#' @export
#' @rdname mpnmf-method
mp_scores <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_scores[index]
}

#' @param n_signatures A scalar integer to specify the number of features to
#' define the program signature.
#' @export
#' @rdname mpnmf-method
mp_signatures <- function(x, s_min = 1 / 3, n_signatures = 20L) {
    mp_scores <- mp_scores(x, s_min = s_min)
    assert_(n_signatures, function(x) {
        is.numeric(x) && is_scalar(x) && x >= 1L
    }, "a positive integer")
    n_signatures <- max(0L, as.integer(n_signatures))
    if (n_signatures > length(mp_scores[[1L]])) {
        cli::cli_warn(paste(
            "{.arg n_signatures} must be smaller than the number of features",
            "in NMF basis matrix"
        ))
        n_signatures <- length(mp_scores[[1L]])
    }
    lapply(mp_scores, function(mp_score) {
        mp_score <- sort(mp_score, decreasing = TRUE)
        names(mp_score)[seq_len(n_signatures)]
    })
}

#' @export
#' @rdname mpnmf-method
mp_samples <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_samples[index]
}

#' @export
#' @rdname mpnmf-method
mp_coverage <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    mp_totals <- lengths(lapply(x$mp_samples, unique))
    mp_coverage <- mp_totals /
        length(unique(unlist(x$mp_samples, FALSE, FALSE)))
    mp_coverage[index]
}

mp_index <- function(x, s_min, call = rlang::caller_call()) {
    assert_s3_class(x, "mpnmf", call = call)
    assert_number(s_min, null_ok = TRUE, call = call)
    mp_samples <- x$mp_samples
    index <- seq_along(mp_samples)
    if (!is.null(s_min)) {
        n_samples <- lengths(lapply(mp_samples, unique))
        if (s_min < 1L && s_min > 0L) {
            n <- length(unique(unlist(mp_samples, FALSE, FALSE)))
            index <- index[(n_samples / n) >= s_min]
        } else {
            index <- index[n_samples >= s_min]
        }
    }
    index
}
