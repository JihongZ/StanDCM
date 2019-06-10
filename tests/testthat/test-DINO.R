library(testthat)
context("test DINO models and model comparsion")

#test_check("StanDCM")


#load("../../data/data.RData")


test_that("run DINO model",{
  Qmatrix2 <- cbind(Qmatrix,rep(1,9))
  Qmatrix2[1,1]<-0
  mod1<-StanDINO.run(Qmatrix2,respMatrix, iter=20,init.list='cdm')
  summary(mod1)
})
#
# test_that("Comparsion between right and wrong model",{
#   mod1 <- StanDINO.run(Qmatrix,respMatrix, iter=20,init.list='cdm')
#   mod2 <- StanDINO.run(misspecifiedQmatrix,respMatrix, iter=20,init.list='cdm')
#   library(loo)
#   loo1 <- loo(mod1)
#   loo2 <- loo(mod2)
#   loo::compare(loo1, loo2)
# })

# test_that("PPMC test",{
#   mod1 <- StanDINO.run(Qmatrix,respMatrix, iter=100,init.list='cdm', warmup = 20)
#   ppmc.dino <- StanDCM.ppmc(mod1, respMatrix, n.sim = 1000, n.burnin = 1, plot.option = TRUE)
#   print(ppmc.dino)
# })


