% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{mp_heatmap}
\alias{mp_heatmap}
\title{Meta program similarity heatmap}
\usage{
mp_heatmap(
  x,
  ...,
  highlight = TRUE,
  s_min = 1/3,
  signatures = NULL,
  textbox_side = NULL,
  textbox_params = list(),
  layer_fun = NULL,
  col = viridisLite::viridis(100, option = "A", direction = -1),
  palette = NULL,
  show_row_names = FALSE,
  show_column_names = FALSE
)
}
\arguments{
\item{x}{A \link[=mp]{mpnmf} object.}

\item{...}{Additional arguments passed to \link[ComplexHeatmap:Heatmap]{Heatmap}.}

\item{highlight}{A boolean value indicates whether highlight significant meta
programs or a \link[grid:gpar]{gpar} object define the highlight
\link[grid:grid.rect]{rect}.}

\item{s_min}{A scalar integer indicates the minimal number of samples or a
scalar numeric (\verb{0 < s_min < 1}) indicates the minimal proportion of samples
to define the meta program. If \code{NULL}, no filters will be applied. Default:
\code{1/3}.}

\item{signatures}{A list of character define the program signatures.}

\item{textbox_side}{A string, indicates the text box side. \code{"left"} or
\code{"right"}. Default: \code{"right"}.}

\item{textbox_params}{Additional arguments passed to
\link[ComplexHeatmap:anno_textbox]{anno_textbox}.}

\item{layer_fun}{Similar as \code{cell_fun}, but is vectorized. Check \url{https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#customize-the-heatmap-body} .}

\item{col}{A vector of colors if the color mapping is discrete or a color mapping  function if the matrix is continuous numbers (should be generated by \code{\link[circlize]{colorRamp2}}). If the matrix is continuous, the value can also be a vector of colors so that colors can be interpolated. Pass to \code{\link[ComplexHeatmap]{ColorMapping}}. For more details and examples, please refer to \url{https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#colors} .}

\item{palette}{A character define to meta program color bar in top and
left. Can also be a string of palette name in
\link[RColorBrewer:ColorBrewer]{RColorBrewer}.}

\item{show_row_names}{Whether show row names.}

\item{show_column_names}{Whether show column names.}
}
\description{
Meta program similarity heatmap
}
