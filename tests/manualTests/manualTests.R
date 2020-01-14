# Tests to be run by hand

# takes about 30 min
allsdy <- CreateConnection("")

try_gei <- function(con) {
  tryCatch(
    con$getGEInputs(),
    warning = function(w) {
      return(w)
    },
    error = function(e) {
      return(e)
    }
  )
}

test_that("returns error if run at project level", {
  res <- try_gei(allsdy)
  expect_true("project" %in% strsplit(res$message, split = " ")[[1]])
})

# Run on server to be local
allsdy$.__enclos_env__$private$.test_files()
