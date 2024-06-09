#' Identify meta programs.
#'
#' @details
#' If `repr` is `tree`, [hclust][stats::hclust] will be used to define the tree
#' structrue, then [cutree][stats::cutree] or
#' [cutreeDynamic][dynamicTreeCut::cutreeDynamic] will be used to cut the tree
#' into groups, in this way, cluster must be a string to define hclust method or
#' a function accepts a [dist][stats::dist] object and return a
#' [hclust][stats::hclust] object. Otherwise, cluster must be the suffix of any
#' igraph community detection algorithm or a function accepts a
#' [graph][igraph::graph] object and return a [communities][igraph::communities]
#' object. For example, cluster="louvain" will use
#' [cluster_louvain][igraph::cluster_louvain].
#'
#' @param nmf_factors A list of NMF factor matrix for each sample, see
#' [basis][NMF::basis].
#' @param n_signatures A scalar integer to specify the number of features to
#' define the program signature.
#' @param ids Sample identifiers, must be the same length of nmf_factors.
#' @param cor_method A character string indicating which correlation coefficient
#' (or covariance) is to be computed. One of "pearson" (default), "kendall", or
#' "spearman": can be abbreviated. See [cor].
#' @param cor_min Program-program similarity were filtered out if their
#' connections were smaller (or equal) than `cor_min` individual tumor modules.
#' Default: `0.3`.
#' @param s_min A scalar integer indicates the minimal number of samples to
#' define the meta program.
#' @param repr A character string of "tree" or "graph" indicates the structure
#' representation of Program-program similarity matrix.
#' @param cluster A character string or function indicating how to clustering
#' the Program-program similarity matrix. Default: if repr is "tree",
#' cluster will be "ward.D2", if repr is "graph", cluster will be "infomap".
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
#' @return A `mpnmf` object.
#' @export
mp <- function(nmf_factors, n_signatures = 20L,
               cor_method = "pearson", cor_min = 0.3, s_min = 3L,
               repr = "tree", cluster = NULL, dynamic = FALSE,
               ..., ids = NULL) {
    assert_(nmf_factors, function(x) {
        is.list(x) && all(vapply(x, is.matrix, logical(1L)))
    }, "a list of NMF factor matrix")

    assert_(n_signatures, function(x) {
        is.numeric(x) && is_scalar(x) && x >= 1L
    }, "a positive integer")
    n_signatures <- as.integer(n_signatures)
    assert_number(cor_min)
    assert_number(s_min)
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
    sample_nms <- rep(sample_nms, times = lengths(nmf_factors))
    program_scores <- lapply(nmf_factors, function(basis) {
        # convert matrix into list ----------------------------
        if (is.null(colnames(basis))) {
            colnames(basis) <- seq_len(ncol(basis))
        }
        apply(basis, 2L, identity, simplify = FALSE)
    })

    # names: {sample}.{nmf} -------------------------------
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
            } else if (!is_scalar(args$h)) {
                cli::cli_abort(
                    "{.arg h} must be a scalar to define the extract tree groups"
                )
            } else if (!is_scalar(args$k)) {
                cli::cli_abort(
                    "{.arg k} must be a scalar to define the extract tree groups"
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
                paste0("cluster_", cluster),
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
                "{.arg cluster} must be a function a character string"
            )
        }
        members <- factor(igraph::membership(comm))
        stats <- list(graph = g, communities = comm)
    }
    # members: a named factor, value is the group identifier, and the name if
    #          the program index
    mp_n_sample <- lengths(lapply(split(sample_nms, members), unique))
    mp_index <- split(names(members), paste0("meta_program_", members))
    mp_index <- lapply(mp_index[mp_n_sample >= s_min], as.integer)

    mp_members <- lapply(mp_index, function(index) program_nms[index])
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
        mp_members,
        similarity = similarity,
        stats = stats,
        mp_scores = mp_scores,
        mp_signatures = lapply(mp_scores, function(mp_score) {
            mp_score <- sort(mp_score, decreasing = TRUE)
            names(mp_score)[seq_len(n_signatures)]
        }),
        class = "mpnmf"
    )
}


#' @export
print.mpnmf <- function(x, ...) {
    for (at in setdiff(names(attributes(x)), "names")) {
        attr(x, at) <- NULL
    }
    print(x)
}
