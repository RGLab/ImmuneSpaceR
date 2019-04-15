#' Save/Load an ImmuneSpaceConnection object from disk
#'
#' Connection can hold a lot of data in cache. If a lot of work has
#' been done (e.g: lots of downloaded datasets and gene-expression matrices),
#' it can be useful to save the connection for later work or even offline use.
#'
#' @param file The file name to be saved to or loaded from
#'
#' @examples
#' # Sample saved connection with pre-downloaded expression matrices and datasets
#' saved <- system.file("extdata/saved_con.rds", package = "ImmuneSpaceR")
#' new_con <- loadConnection(saved)
#' new_con
#' names(new_con$cache)
#' \dontrun{
#' saveConnection(new_con, tempfile())
#' }
#' 
#' @rdname loadConnection
#' @export
#' @return An ImmuneSpaceConnection object
loadConnection <- function(file) {
  con <- readRDS(file = file)
  conType <- class(con)

  if (conType == "ImmuneSpaceConnection") {
    labkey.url.base <- con$config$labkey.url.base
  } else {
    stop("invalid ImmuneSpaceConnection object!")
  }

  # init labkey.setCurlOptions
  labkey.setCurlOptions(ssl_verifyhost = 2, sslversion = 1)
  con
}


#' @param con An \code{ImmuneSpaceConnection}. The connection to save to file.
#'  To be loaded later using \code{loadConnection}.
#'
#' @rdname loadConnection
#' @export
saveConnection <- function(con, file) {
  saveRDS(con, file = file)
}
