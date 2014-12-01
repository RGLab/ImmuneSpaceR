context("test suit for multi-studies connection (SDY269, SDY180, SDY61)")

test_that("studyNames",{
      expect_equal(multiSdy$studyNames(), studies)
    })

test_that("listDatasets",{
      thisRes <- capture.output(multiSdy$listDatasets())
      expect_equal(expectRes[["listDS_multi"]], thisRes)
      
    })

test_that("listGEAnalysis",{
      thisRes <- multiSdy$listGEAnalysis()
      name <- "listGE_multi"
      expect_equal(expectRes[[name]], thisRes)
    })


test_that("getDataset",{
      
      thisRes <- multiSdy$getDataset("hai")
      name <- "hai_multi"
      expect_equivalent(expectRes[[name]], thisRes)
      
      thisRes <- multiSdy$getDataset("elispot")
      name <- "elispot_multi"
      expect_equivalent(expectRes[[name]], thisRes)
      
      
    })

test_that("getGEMatrix",{
      thisRes <- multiSdy$getGEMatrix(c("TIV_2007", "Saline_group1", "LAIV_2008"), summary = TRUE)
      name <- "getGE_multi"
      expect_equal(expectRes[[name]], thisRes)
            
    })



