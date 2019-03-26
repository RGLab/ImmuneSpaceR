#' template_IS
#'
#' A HTML template for knitted reports that matches ImmuneSpace's graphic style.
#' It is based on \code{\link[rmarkdown]{html_document}} from the \pkg{rmarkdown}
#' package with css, theme, and template parameters disabled.
#'
#' @param ... See \code{\link[rmarkdown]{html_document}}
#'
#' @return R Markdown output format to pass to \code{\link[rmarkdown]{render}}
#'
#' @details
#' See the documentation for \code{\link[rmarkdown]{html_document}} or the
#' \href{http://rmarkdown.rstudio.com/html_document_format.html}{oneline documentation}
#' for additional details on using the html_document format.
#' Compared to html_document, it:
#' \itemize{
#' \item uses a custom css stylesheet
#' \item does not use bootstrap themes
#' }
#'
#' @examples
#' \dontrun{
#' library(ImmuneSpaceR)
#' rmarkdown::render("input.Rmd", template_IS())
#' rmarkdown::render("input.Rmd", template_IS(toc = TRUE))
#' }
#' template_IS()
#' @importFrom rmarkdown html_document
#' @export
template_IS <- function(...) {
  html_document(
    css = system.file(
      "rmarkdown/templates/ImmuneSpace/resources/IStemplate.css",
      package = "ImmuneSpaceR"
    ),
    theme = NULL,
    template = "default",
    ...
  )
}
