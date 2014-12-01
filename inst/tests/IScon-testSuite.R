context("SDY269 test suite...")

test_that("listDatasets",{
      thisRes <- capture.output(sdy269$listDatasets())
      name <- "listDS_269"
      expect_equal(expectRes[[name]], thisRes)
      
    })




test_that("listGEAnalysis",{
      thisRes <- sdy269$listGEAnalysis()
      name <- "listGE_269"
      expect_equal(expectRes[[name]], thisRes)
      
    })

test_that("getDataset",{
      
      
      thisRes <- sdy269$getDataset("hai")
      thisRes[, study:="SDY269"]
      expect_equal(expectRes[["hai_multi"]][study == "SDY269"], thisRes)
      
      thisRes <- sdy269$getDataset("elispot")
      thisRes[, study:="SDY269"]
      expect_equal(expectRes[["elispot_multi"]][study == "SDY269"], thisRes)
      
      thisRes <- sdy269$getDataset("pcr")
      name <- "pcr_269"
      expect_equal(expectRes[[name]], thisRes)
      
      
      thisRes <- sdy269$getDataset("elisa")
      name <- "elisa_269"
      expect_equal(expectRes[[name]], thisRes)
      
      thisRes <- sdy269$getDataset("fcs_analyzed_result")
      name <- "fcs_analyzed_result_269"
      expect_equal(expectRes[[name]], thisRes)
      
      thisRes <- sdy269$getDataset("demographics")
      name <- "demographics_269"
      expect_equal(expectRes[[name]], thisRes)
            
    })

test_that("getGEMatrix",{
    thisRes <- sdy269$getGEMatrix("LAIV_2008")
    
    if(isSubset)#result of multi-study version is a subset of original
    {
      name <- "getGE_269_subset"
      expect_equal(expectRes[[name]], thisRes)
    }else{
      name <- "getGE_269"
      expect_equal(expectRes[[name]], thisRes)
    }
     
    })

## we need to either make it return ggplot object or separate the main logic from the quick_plot
## before writing test for quick_plot
#test_that("quick_plot",{
#      expect_equal('', digest(sdy269$quick_plot("hai")))
#    })


