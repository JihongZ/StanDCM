library(testthat)
context("test LCDM and DINA")

#test_check("StanDCM")


#load("../../data/data.RData")

# test_that("runing LCDM model", {
#   StanLCDM.run(Qmatrix,respMatrix,iter=20,init.list='cdm',chains = 3)
# })

test_that("test Loo fit model comparision", {
  estimatedMod<-StanLCDM.run(Qmatrix,respMatrix,iter=20,init.list='cdm',chains = 3)
  StanDCM::StanLCDM.loofit(estimatedMod)
})


# test_that("run DINA model",{
#   Qmatrix2 <- cbind(Qmatrix,rep(1,9))
#   Qmatrix2[1,1]<-0
#   StanDINA.run(Qmatrix2,respMatrix,iter=20,init.list='cdm')
# })
