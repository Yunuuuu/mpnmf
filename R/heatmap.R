#' Meta program similarity heatmap
#'
#' @param x A [mpnmf][mp] object.
#' @param ... Additional arguments passed to [Heatmap][ComplexHeatmap::Heatmap].
#' @param highlight A boolean value indicates whether highlight significant meta
#' programs or a [gpar][grid::gpar] object define the highlight
#' [rect][grid::grid.rect].
#' @inheritParams mpnmf-method
#' @inheritParams ComplexHeatmap::Heatmap
#' @param palette A named character define to meta program color bar in top and
#' left. Can also be a string of palette name in
#' [RColorBrewer][RColorBrewer::brewer.pal].
#' @export
mp_heatmap <- function(
    x, ..., highlight = TRUE, s_min = 1 / 3, layer_fun = NULL,
    col = viridisLite::viridis(100, option = "A", direction = -1),
    palette = NULL, show_row_names = FALSE, show_column_names = FALSE) {
    assert_s3_class(x, "mpnmf")
    if (isTRUE(highlight)) {
        highlight <- grid::gpar(lty = 2, lwd = 2, col = "red")
    } else if (inherits(highlight, "gpar")) {
        highlight$fill <- NA
    } else if (isFALSE(highlight)) {
        cli::cli_abort(
            "{.arg highlight} must be a boolean value or a {.cls gpar} object"
        )
    }
    mat <- attr(x, "similarity")
    stats <- attr(x, "stats")
    n <- nlevels(stats$members)
    if (is.null(palette)) {
        if (n > 12L) {
            palette <- ComplexHeatmap:::default_col(stats$members)
        } else {
            palette <- RColorBrewer::brewer.pal(n, "Paired")
            names(palette) <- levels(stats$members)
        }
    } else if (rlang::is_string(palette) && is.null(names(palette))) {
        palette <- match.arg(palette, rownames(RColorBrewer::brewer.pal.info))
        max <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
        if (max < n) {
            cli::cli_abort(
                "{palette} is not sufficient to provide {.val {n}} colors"
            )
        }
        palette <- RColorBrewer::brewer.pal(n, palette)
        names(palette) <- levels(stats$members)
    }
    if (inherits(highlight, "gpar")) {
        sig_programs <- mp_programs(x, s_min = s_min)
        sig_index <- match(levels(stats$members), names(sig_programs))
        highlight_layer_fun <- function(highlight) {
            slice_numbers <- strsplit(grid::current.viewport()$name,
                "_",
                fixed = TRUE
            )[[1L]]
            slice_numbers <- as.integer(utils::tail(slice_numbers, n = 2L))
            if (any(.subset(slice_numbers, 1L) == sig_index) &&
                # the second and third items are the slice row and column number
                .subset(slice_numbers, 1L) == .subset(slice_numbers, 2L)) {
                grid::grid.rect(gp = highlight)
            }
        }
        if (is.null(layer_fun)) {
            internal_layer_fun <- function(j, i, x, y, w, h, fill) {
                highlight_layer_fun(highlight)
            }
        } else {
            internal_layer_fun <- function(j, i, x, y, w, h, fill) {
                layer_fun(j, i, x, y, w, h, fill)
                highlight_layer_fun(highlight)
            }
        }
    } else {
        internal_layer_fun <- layer_fun
    }
    dist_fn <- function(x) stats::as.dist(1 - x[, rownames(x)])

    args <- list(
        matrix = mat,
        left_annotation = ComplexHeatmap::HeatmapAnnotation(
            MP = stats$members, col = list(MP = palette),
            which = "row", show_legend = FALSE,
            show_annotation_name = FALSE
        ),
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
            MP = stats$members, col = list(MP = palette),
            which = "column",
            show_annotation_name = FALSE,
            annotation_legend_param = list(title = NULL)
        ),
        col = col, ...,
        layer_fun = internal_layer_fun,
        row_split = stats$members,
        column_split = stats$members,
        clustering_distance_rows = dist_fn,
        clustering_distance_columns = dist_fn,
        show_row_names = show_row_names,
        show_column_names = show_column_names
    )
    args$name <- args$name %||% "Similarity"
    if (!rlang::has_name(args, "row_title")) {
        args["row_title"] <- list(NULL)
    }
    if (!rlang::has_name(args, "column_title")) {
        args["column_title"] <- list(NULL)
    }
    if (!is.null(stats$tree) && !is.null(stats$tree$method)) {
        # restore tree object from `mpnmf` object
        args$clustering_method_rows <- args$clustering_method_columns <-
            stats$tree$method
    } else {
        args$clustering_method_rows <- args$clustering_method_rows %||%
            "ward.D2"
        args$clustering_method_columns <- args$clustering_method_columns %||%
            "ward.D2"
    }
    do.call(ComplexHeatmap::Heatmap, args)
}
