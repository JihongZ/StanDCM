library(testthat)
context("Test multigroup DINO")

test_that("multigoup DINO run",{
  # save.name="DINO_uninf_multiG"
  # Qmatrix=Qmatrix;response.matrix=respMatrix;GroupID=c(rep(1,500),rep(2,500));
  # fixeditem.vector=c(1,3);class.equal=F;chain.num=3;init.list='random';
  # save.path=getwd();
  # warmup=10;control.list<-list(adapt_delta=0.82);iter=200
  options(mc.cores = parallel::detectCores())
  StanDINO_mG.script(Qmatrix,2,c(1,2),class.equal=F)
  StanDINO_mG.run(Qmatrix = Qmatrix,
                  response.matrix= respMatrix,
                  GroupID = c(rep(1,500),rep(2,500)),
                  fixeditem.vector = c(1,3),
                  class.equal = F,
                  script.path=NA,
                  save.path=getwd(),
                  save.name = "DINO_uninf_multiG",
                  iter = 200,
                  warmup = 10,
                  chain.num = 3,
                  init.list = 'random',
                  control.list = list(adapt_delta=0.82))

  })
