library(testthat)
context("test LCDM and DINA")

#test_check("StanDCM")


#load("../../data/data.RData")


test_that("run DINA model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  StanDINA.run(Qmatrix2, respMatrix, iter=20, init.list='cdm')
})


