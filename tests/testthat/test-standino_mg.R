library(testthat)
context("Test multigroup DINO")


test_that("mg run",{
  GroupID = c(rep(1, 500), rep(2, 500))
  StanDINO_mG.run(Qmatrix = Qmatrix, response.matrix = respMatrix, GroupID = GroupID, fixeditem.vector = c(2,3),
                  class.equal = T)
})
