## ----echo=TRUE-----------------------------------------------------------
labkey.url.base="https://www.immunsepace.org"
labkey.url.path="Studies/SDY269"
labkey.user.email="unknown_user at not_a_real_domain.org"

## ----,message=FALSE------------------------------------------------------
require(ImmuneSpaceR)
labkey.url.base="https://www.immunespace.org" 
labkey.url.path="Studies/SDY269"
labkey.email.user="gfinak at fhcrc.org"

## ----, CreateConnection--------------------------------------------------
sdy269<-CreateConnection(study="SDY269")
sdy269

## ----getDataset----------------------------------------------------------
sdy269$getDataset("hai")

## ----getGEMatrix---------------------------------------------------------
sdy269$getGEMatrix("LAIV_2008")

## ----getGEMatrix-multiple------------------------------------------------
sdy269$getGEMatrix(c("TIV_2008", "LAIV_2008"))

## ------------------------------------------------------------------------
result <- sdy269$getDataset("mickey_mouse")

## ----, dev='CairoPNG'----------------------------------------------------
sdy269$quick_plot("hai")

sdy269$quick_plot("elisa")

