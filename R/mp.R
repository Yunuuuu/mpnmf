#' Identify meta programs.
#'
#' @details
#' If `repr` is `tree`, [hclust][stats::hclust] will be used to define the tree
#' structrue, then [cutree][stats::cutree] or
#' [cutreeDynamic][dynamicTreeCut::cutreeDynamic] will be used to cut the tree
#' into groups, in this way, `cluster` must be a string to define hclust method
#' or a function accepts a [dist][stats::dist] object and return a
#' [hclust][stats::hclust] object. Otherwise, `cluster` must be the suffix of
#' any igraph community detection algorithm or a function accepts a
#' [graph][igraph::graph] object and return a [communities][igraph::communities]
#' object. For example, cluster="louvain" will use
#' [cluster_louvain][igraph::cluster_louvain].
#'
#' @param nmf_factors A list of NMF factor matrix for each sample, see
#' [basis][NMF::basis].
#' @param ids Sample identifiers, must be the same length of nmf_factors.
#' @param cor_method A character string indicating which correlation coefficient
#' (or covariance) is to be computed. One of "pearson" (default), "kendall", or
#' "spearman": can be abbreviated. Details see [cor].
#' @param cor_min Program-program similarity were filtered out if their
#' connections were smaller (or equal) than `cor_min` individual tumor modules.
#' Default: `0.3`.
#' @param repr A character string of "tree" or "graph" indicates the structure
#' representation of Program-program similarity matrix.
#' @param cluster A character string or function indicating how to clustering
#' the program-program similarity matrix. Default: if repr is `"tree"`,
#' cluster will be `"ward.D2"`, if repr is `"graph"`, cluster will be
#' `"infomap"`.
#' @param dynamic A boolean value indicates whether to use
#' [cutreeDynamic][dynamicTreeCut::cutreeDynamic] to define the tree groups.
#' Only be used when repr is "tree".
#' @param ... Additional arguments passed to [cutree][stats::cutree] or
#' [cutreeDynamic][dynamicTreeCut::cutreeDynamic] or `igraph::cluster_*`
#' function (when `cluster` was passed as a function, ...  will not be used).
#' @note When using a cutom function in `cluster`, you must follow the tree
#' nodes (or graph vertex) names, that means you must return the groups in the
#' same order of the tree nodes (or graph vertex) name. Since the internal will
#' restore the program names using the tree nodes (or graph vertex) name.
#' @return
#' `mp`: A list of class `mpnmf` object with following elements:
#'
#'  - `mp_programs`: A list of meta programs.
#'  - `mp_samples`: A list of samples define the meta programs.
#'  - `mp_scores`: A list of meta program scores, which were defined as the mean
#'                 NMF factor basis across the component programs.
#'
#' In addition, following attributes are attached with this object
#'  - `similarity`: similarity matrix.
#'  - `stats`: A list of statistics for tree or graph object
#'
#' `mp_signatures`: A list of features to define the meta program signature.
#' `mp_samples`: A list of samples indicates where the meta program is.
#' `mp_coverage`: A numeric indicates the fraction of samples the meta program
#'                was detected
#'
#' @export
mp <- function(nmf_factors,
               cor_method = "pearson", cor_min = 0.3,
               repr = "tree", cluster = NULL, dynamic = FALSE,
               ..., ids = NULL) {
    assert_(nmf_factors, function(x) {
        is.list(x) && all(vapply(x, is.matrix, logical(1L)))
    }, "a list of NMF factor matrix")

    assert_number(cor_min)
    repr <- match.arg(repr, c("tree", "graph"))

    # check feature name provided ------------------------
    all_features <- lapply(nmf_factors, rownames)
    if (any(vapply(all_features, is.null, logical(1L)))) {
        cli::cli_abort(
            "all basis matrix must have rownames to define the features"
        )
    }
    features <- all_features[[1L]]
    for (i in seq_along(all_features)) {
        if (!identical(features, all_features[[i]])) {
            cli::cli_abort(
                "all basis matrix must have identical rownames (features)"
            )
        }
    }

    # add sample names if no name provided  ------------
    if (!is.null(ids)) {
        if (length(ids) != length(nmf_factors)) {
            cli::cli_abort(
                "{.arg id} must be character with the same length of {.arg nmf_factors}"
            )
        }
        names(nmf_factors) <- ids
    }

    sample_nms <- names(nmf_factors)
    if (is.null(sample_nms)) {
        names(nmf_factors) <- seq_along(nmf_factors)
    } else if (anyDuplicated(sample_nms)) {
        cli::cli_abort("names of {.arg nmf_factors} must be unique")
    } else if (anyNA(sample_nms) || any(sample_nms == "")) {
        cli::cli_abort("names of {.arg nmf_factors} cannot be missing")
    }

    # convert NMF basis matrix into list -----------------------
    program_scores <- lapply(nmf_factors, function(basis) {
        # convert matrix into list ----------------------------
        if (is.null(colnames(basis))) {
            colnames(basis) <- seq_len(ncol(basis))
        }
        apply(basis, 2L, identity, simplify = FALSE)
    })

    # flatten vectors -------------------------------
    sample_nms <- rep(sample_nms, times = lengths(program_scores))
    # names: {sample}.{nmf}
    program_scores <- unlist(program_scores, recursive = FALSE)
    program_nms <- names(program_scores)
    duplicates <- duplicated(program_nms) |
        duplicated(program_nms, fromLast = TRUE)
    if (any(duplicates)) {
        cli::cli_warn(paste(
            "Found duplicated program names, this should be due to",
            "the exist of duplicated colnames in the NMF factor basis",
            "of sample{?s} ({unique(sample_nms[duplicates])})"
        ))
    }
    program_index <- seq_along(program_scores)

    # we use index as the program names in case of duplicated names
    names(program_scores) <- program_index
    similarity <- stats::cor(
        do.call(base::cbind, program_scores),
        method = cor_method
    )
    similarity[abs(similarity) <= cor_min] <- 0

    # cluster all programs to identify meta programs -----------
    if (repr == "tree") {
        assert_bool(dynamic)
        dist <- stats::as.dist(1 - similarity)
        cluster <- cluster %||% "ward.D2"
        if (is.character(cluster)) {
            if (!is_scalar(cluster)) {
                cli::cli_abort("{.arg cluster} must be a scalar string")
            }
            hcl <- do.call(stats::hclust, list(d = dist, method = cluster))
        } else if (is.function(cluster)) {
            if (!inherits(hcl <- cluster(dist), "hclust")) {
                cli::cli_abort(
                    "{.arg cluster} must return a {.cls hclust} object"
                )
            }
        } else {
            cli::cli_abort(
                paste(
                    "{.arg cluster} must be a string to define the `hclust`",
                    "method or a function accepts a `dist` object",
                    "and returns a `hclust` tree"
                )
            )
        }
        stats <- list(tree = hcl)
        if (dynamic) {
            members <- dynamicTreeCut::cutreeDynamic(
                dendro = hcl, ...,
                distM = as.matrix(dist)
            )
        } else {
            args <- list(...)
            if (is.null(args$k) && is.null(args$h)) {
                args$h <- max(hcl$height) / 2
            } else if (!is.null(args$h) && !is_scalar(args$h)) {
                cli::cli_abort(
                    "{.arg h} must be a scalar to define the exact tree groups"
                )
            } else if (!is.null(args$k) && !is_scalar(args$k)) {
                cli::cli_abort(
                    "{.arg k} must be a scalar to define the exact tree groups"
                )
            }
            members <- do.call(stats::cutree, c(list(tree = hcl), args))
        }
        members <- factor(members)
    } else {
        # vertex names will be the same with matrix names
        g <- igraph::graph_from_adjacency_matrix(similarity,
            mode = "undirected", diag = FALSE, weighted = TRUE
        )
        cluster <- cluster %||% "infomap"
        if (is.character(cluster)) {
            if (!is_scalar(cluster)) {
                cli::cli_abort("{.arg cluster} must be a scalar string")
            }
            cluster <- utils::getFromNamespace(
                paste("cluster", cluster, sep = "_"),
                ns = "igraph"
            )
            comm <- cluster(graph = g, ...)
        } else if (is.function(cluster)) {
            if (!inherits(comm <- cluster(graph = g), "communities")) {
                cli::cli_abort(
                    "{.arg cluster} must return a {.cls communities} object"
                )
            }
        } else {
            cli::cli_abort(
                paste(
                    "{.arg cluster} must be a string to define the `cluster_*`",
                    "method or a function accepts a `graph` object",
                    "and returns a `communities` object"
                )
            )
        }
        members <- factor(igraph::membership(comm))
        stats <- list(graph = g, communities = comm)
    }

    # members: a named factor, value is the group identifier, and the name is
    #          the program index
    levels(members) <- paste0("MP", levels(members))
    mp_index <- lapply(split(names(members), members), as.integer)

    mp_programs <- lapply(mp_index, function(index) program_nms[index])
    mp_samples <- lapply(mp_index, function(index) sample_nms[index])
    mp_scores <- lapply(mp_index, function(index) {
        component_scores <- program_scores[index]
        mp_score <- rowMeans(do.call(base::cbind, component_scores))
        mp_score / sum(mp_score)
    })

    # restore names -----------------------------------------
    names(members) <- program_nms
    dimnames(similarity) <- list(program_nms, program_nms)
    stats <- c(stats, list(members = members))

    # return results -----------------------------------------
    structure(
        list(
            mp_programs = mp_programs,
            mp_samples = mp_samples,
            mp_scores = mp_scores
        ),
        similarity = similarity,
        stats = stats,
        class = "mpnmf"
    )
}

#' @param x A `mpnmf` object.
#' @param s_min A scalar integer indicates the minimal number of samples or a
#' scalar numeric (`0 < s_min < 1`) indicates the minimal proportion of samples
#' to define the meta program. If `NULL`, no filters will be applied. Default:
#' `1/3`.
#' @export
#' @rdname mp
print.mpnmf <- function(x, s_min = 1 / 3, ...) {
    print(mp_programs(x, s_min))
}

#' @export
#' @rdname mp
mp_programs <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_programs[index]
}

#' @export
#' @rdname mp
mp_scores <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_scores[index]
}

#' @param n_signatures A scalar integer to specify the number of features to
#' define the program signature.
#' @export
#' @rdname mp
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
#' @rdname mp
mp_samples <- function(x, s_min = 1 / 3) {
    index <- mp_index(x, s_min)
    x$mp_samples[index]
}

#' @export
#' @rdname mp
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
