context("test suit for multi-studies connection (SDY269, SDY180, SDY61)")

test_that("studyNames",{
      expect_equal(multiSdy$studyNames(), studies)
    })

test_that("listDatasets",{
      expect_equal('f37d99cebbcc6686e49feebdc8c5af3e', digest(capture.output(multiSdy$listDatasets())))
      
    })

test_that("listGEAnalysis",{
      expect_equal('45bcd4ad93f7d99b2761924b9d15f709', digest(multiSdy$listGEAnalysis()))
    })


test_that("getDataset",{
      expect_equal('b04021a57520c8f3c8246a5b4c28fd2c', digest(multiSdy$getDataset("hai")))
      expect_equal('9c432825a7cbbb59362f63fd480d6e5c', digest(multiSdy$getDataset("elispot")))
      
    })

test_that("getGEMatrix",{
      expect_equal('dc04655832df53d0d2e55c84db404691', digest(multiSdy$getGEMatrix(c("TIV_2007", "Saline_group1", "LAIV_2008"), summary = TRUE)))
    })


