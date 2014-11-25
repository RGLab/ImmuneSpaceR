context("SDY269 test suite...")

test_that("listDatasets",{
      expect_equal('4c741bbb2d5c4863ed0b569758c7d72d', digest(capture.output(sdy269$listDatasets())))
    })



test_that("listGEAnalysis",{
      expect_equal('490f50d315badb143dca6e0c518d27c2', digest(sdy269$listGEAnalysis()))
    })

test_that("getDataset",{
      expect_equal('cd9eb24ee29221a438a9be2f80fb14ff', digest(sdy269$getDataset("hai")))
      expect_equal('acce67e4d065e5066f2227f6082f2d7d', digest(sdy269$getDataset("elispot")))
      expect_equal('1193139c6221fb3b6251a64c9736227a', digest(sdy269$getDataset("pcr")))
      expect_equal('cd4224d26715c50007f0f418a6fa022e', digest(sdy269$getDataset("elisa")))
      expect_equal('5fb99008a8412a43200eae4f82c9cd4e', digest(sdy269$getDataset("fcs_analyzed_result")))
      expect_equal('f6722b1efaed6887778e2d3182b7c300', digest(sdy269$getDataset("demographics")))
      
    })

test_that("getGEMatrix",{
    thisRes <- sdy269$getGEMatrix("LAIV_2008")
    
    if(isSubset)#result of multi-study version is a subset of original 
      expect_equal('427f18ce61c15946b906769f823a1f23', digest(thisRes))
    else
      expect_equal('36beb0adde3e48af1f93aa4916de2e0e', digest(thisRes))
      
    })

## we need to either make it return ggplot object or separate the main logic from the quick_plot
## before writing test for quick_plot
#test_that("quick_plot",{
#      expect_equal('', digest(sdy269$quick_plot("hai")))
#    })


