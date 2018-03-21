#' theme_IS
#'
#' Theme that matches ImmuneSpace's graphic style. The theme modifies the
#' background, the grid lines, the axis, and the colors used by continuous and
#' gradient scales.
#'
#' @param base_size A \code{numeric}. Base font size.
#'
#' @return A theme object
#'
#' @details
#' List of modified ggplot2 elements: \code{panel.background},
#' \code{panel.grid.major}, \code{panel.grid.minor}, \code{axis.ticks},
#' \code{axis.line.x}, \code{axis.line.y}, \code{plot.title}, and
#' \code{strip.background}.
#'
#' The default \code{scale_fill_gradient}, \code{scale_fill_continuous},
#' \code{scale_colour_gradient} and \code{scale_colour_continous} are also
#' replaced by a custom scale.
#'
#' @importFrom ggplot2 theme theme_classic element_line element_text element_rect
#' @importFrom ggplot2 continuous_scale rel
#' @export
#' @examples
#' library(ggplot2)
#' p <- ggplot(data = mtcars) + geom_point(aes(x = mpg, y = cyl, color = hp)) + facet_grid(vs ~ am)
#' p + theme_IS()
theme_IS <- function(base_size = 12) {
  .override_scale()
  theme(text = element_text(size = base_size)) +
    #theme(axis.text.x = element_text(angle = 45)) +
    theme(panel.background = element_rect(fill = "#FDF6E3")) +
    theme(panel.grid.major = element_line(colour = "#ded8d5", linetype = "dashed")) +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.ticks = element_blank()) +
    #theme(axis.line = element_line(size = 0.5, colour = "black")) +
    theme(axis.line.x = element_line(size = 0.5, colour = "black")) +
    theme(axis.line.y = element_line(size = 0.5, colour = "black")) +
    theme(plot.title = element_text(size = rel(1))) +
    theme(strip.background = element_rect(colour = "white", fill = "white")) #bg of facets
}


#' ImmuneSpace palette
#'
#' Create a color gradient of the selected length that matches the ImmuneSpace
#' theme.
#'
#' @param n A \code{numeric}. The length of the desired palette.
#' @return A \code{character} vector colors in hexadecimal code of length
#'  \code{n}.
#'
#' @importFrom gplots colorpanel
#' @export
#' @examples
#' plot(1:10, col = ISpalette(10), cex = 10, pch = 16)
ISpalette <- function(n) {
  colorpanel(n, low = "#268bd2", mid = "#fdf6e3", high = "#dc322f")
}


# Add scales that use a different color scheme to the environment.
#' @importFrom ggplot2 update_geom_defaults
#' @importFrom scales seq_gradient_pal
.override_scale <- function(envir = as.environment(1)) {
  update_geom_defaults("boxplot", list(fill = "#268bd2"))
  scale_updates <- list(
    scale_fill_continuous = function(...) continuous_scale('fill', 'scale_IS', seq_gradient_pal("#268bd2", "#dc322f"), ...),
    scale_fill_gradient = function(...) continuous_scale('fill', 'scale_IS', seq_gradient_pal("#268bd2", "#dc322f"), ...),
    scale_colour_gradient = function(...) continuous_scale('colour', 'scale_IS', seq_gradient_pal("#268bd2", "#dc322f"), ...),
    scale_colour_continuous = function(...) continuous_scale('colour', 'scale_IS', seq_gradient_pal("#268bd2", "#dc322f"), ...)
  )

  Map(
    function(name, f) assign(name, f, envir = envir),
    names(scale_updates),
    scale_updates
  )
}
