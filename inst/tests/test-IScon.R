context("ImmuneSpace Connection")

#' create multi studies from scratch
studies <- c("SDY269", "SDY180", "SDY61")

test_that("CreateConnection -- multiple studies",{
      
      multiSdy <<- CreateConnection(study = studies)
      expect_is(multiSdy, "ImmuneSpaceConnectionList")
    })
#' test API of con list
source("ISconList-testSuite.R", local = TRUE)

#' save/load of multi study to/from disck
test_that("save/load ConnectionList",{
      tmp <- tempfile()
      saveConnection(multiSdy, file = tmp)
      multiSdy <<- loadConnection(tmp)
      expect_is(multiSdy, "ImmuneSpaceConnectionList")
      
    })
source("ISconList-testSuite.R", local = TRUE)

#' single study
test_that("subset single study",{
  sdy269 <<- multiSdy$study("SDY269") 
  expect_is(sdy269, "ImmuneSpaceConnection")
})

isSubset <- TRUE
source("IScon-testSuite.R", local = TRUE)

#' save/load of single study to/from disck
test_that("save/load Connections",{
      tmp <- tempfile()
      saveConnection(sdy269, file = tmp)
      sdy269 <- loadConnection(tmp)
      expect_is(sdy269, "ImmuneSpaceConnection")
      
    })
isSubset <- TRUE
source("IScon-testSuite.R", local = TRUE)

#' create single study from scratch
test_that("CreateConnection",{
      sdy269 <<- CreateConnection(study= c("SDY269"))
      expect_is(sdy269, "ImmuneSpaceConnection")
    })
isSubset <- FALSE
source("IScon-testSuite.R", local = TRUE)
