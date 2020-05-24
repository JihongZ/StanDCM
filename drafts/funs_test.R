library(StanDCM)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


####################### Use default Qmatrix and Data ##############################
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 4)
summary(summary(mod.LCDM)$summary[,"Rhat"])

mod.CRUM <- StanCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(summary(mod.CRUM)$summary[,"Rhat"])

mod.DINA <- StanDINA.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(summary(mod.DINA)$summary[,"Rhat"])

mod.DINO <- StanDINO.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(summary(mod.DINO)$summary[,"Rhat"])

mod.NCRUM <- StanNCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(summary(mod.NCRUM)$summary[,"Rhat"])

mod.ORDM <- StanNCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
summary(summary(mod.ORDM)$summary[,"Rhat"])


####################### Use wrong Qmatrix and original Data ##############################
Qmatrix = cbind(Qmatrix, c(1,1,1,1,0,0,0,0,0))

path = "~/Desktop/Out/Condition2/"
#modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
sink(paste0(path, "output.txt"), type = c("output", "message"))
print(Qmatrix)
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.LCDM)$summary[,"Rhat"])


mod.CRUM <- StanCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.CRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.CRUM)$summary[,"Rhat"])

mod.DINA <- StanDINA.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.DINA@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINA)$summary[,"Rhat"])

mod.DINO <- StanDINO.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.DINO@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINO)$summary[,"Rhat"])

mod.NCRUM <- StanNCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.NCRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.NCRUM)$summary[,"Rhat"])

mod.ORDM <- StanORDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 5000, warmup = 1000, chain.num = 5)
modelname =  strsplit(mod.ORDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.ORDM)$summary[,"Rhat"])
sink()


####################### Condition 3 Use Qmatrix and original Data ##############################
set.seed(1222)
a1 = sample(c(0,1), 9, replace = T)
Qmatrix = cbind(Qmatrix, a1)

path = "~/Desktop/Out/Condition2/"
#modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
sink(paste0(path, "output2.LCDM.txt"), type = c("output", "message"))
print(Qmatrix)
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.LCDM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output2.CRUM.txt"), type = c("output", "message"))
mod.CRUM <- StanCRUM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.CRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.CRUM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output2.DINA.txt"), type = c("output", "message"))
mod.DINA <- StanDINA.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINA@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINA)$summary[,"Rhat"])
sink()

sink(paste0(path, "output2.DINO.txt"), type = c("output", "message"))
mod.DINO <- StanDINO.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINO@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINO)$summary[,"Rhat"])
sink()

sink(paste0(path, "output2.NCRUM.txt"), type = c("output", "message"))
mod.NCRUM <- StanNCRUM.run(Qmatrix = Qmatrix,
                           response.matrix = respMatrix,
                           iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.NCRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.NCRUM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output2.ORDM.txt"), type = c("output", "message"))
mod.ORDM <- StanORDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.ORDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.ORDM)$summary[,"Rhat"])
sink()

system("shutdown -h now")

####################### Condition 4 Use Qmatrix and original Data ##############################
Qmatrix4 = Qmatrix
Qmatrix4[5,2] = 1
Qmatrix4[6,2] = 1
Qmatrix4 = Qmatrix4[,1:2]

path = "~/Desktop/Out/Condition4/"
#modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
sink(paste0(path, "output.LCDM.txt"), type = c("output", "message"))
print(Qmatrix)
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix4,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.LCDM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output.CRUM.txt"), type = c("output", "message"))
mod.CRUM <- StanCRUM.run(Qmatrix = Qmatrix4,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.CRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.CRUM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output.DINA.txt"), type = c("output", "message"))
mod.DINA <- StanDINA.run(Qmatrix = Qmatrix4,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINA@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINA)$summary[,"Rhat"])
sink()

sink(paste0(path, "output.DINO.txt"), type = c("output", "message"))
mod.DINO <- StanDINO.run(Qmatrix = Qmatrix4,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINO@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINO)$summary[,"Rhat"])
sink()

sink(paste0(path, "output.NCRUM.txt"), type = c("output", "message"))
mod.NCRUM <- StanNCRUM.run(Qmatrix = Qmatrix4,
                           response.matrix = respMatrix,
                           iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.NCRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.NCRUM)$summary[,"Rhat"])
sink()

sink(paste0(path, "output.ORDM.txt"), type = c("output", "message"))
mod.ORDM <- StanORDM.run(Qmatrix = Qmatrix4,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.ORDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.ORDM)$summary[,"Rhat"])
sink()

system("shutdown -h now")


####################### Condition 5 Use Qmatrix and original Data ##############################
Qmatrix5 = Qmatrix
Qmatrix5[,3] = 1

path = "~/Desktop/Out/Condition5/"
#modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
sink(paste0(path, "output.LCDM.txt"), type = c("output", "message"))
print(Qmatrix)
mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix5,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.LCDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.LCDM)$summary[,"Rhat"])

sink(paste0(path, "output.CRUM.txt"), type = c("output", "message"))
mod.CRUM <- StanCRUM.run(Qmatrix = Qmatrix5,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.CRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.CRUM)$summary[,"Rhat"])

sink(paste0(path, "output.DINA.txt"), type = c("output", "message"))
mod.DINA <- StanDINA.run(Qmatrix = Qmatrix5,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINA@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINA)$summary[,"Rhat"])

sink(paste0(path, "output.DINO.txt"), type = c("output", "message"))
mod.DINO <- StanDINO.run(Qmatrix = Qmatrix5,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.DINO@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.DINO)$summary[,"Rhat"])

sink(paste0(path, "output.NCRUM.txt"), type = c("output", "message"))
mod.NCRUM <- StanNCRUM.run(Qmatrix = Qmatrix5,
                           response.matrix = respMatrix,
                           iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.NCRUM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.NCRUM)$summary[,"Rhat"])

sink(paste0(path, "output.ORDM.txt"), type = c("output", "message"))
mod.ORDM <- StanORDM.run(Qmatrix = Qmatrix5,
                         response.matrix = respMatrix,
                         iter = 3000, warmup = 1000, chain.num = 4)
modelname =  strsplit(mod.ORDM@stanmodel@model_cpp$model_cppname, split="_")[[1]][2]
print(modelname)
summary(summary(mod.ORDM)$summary[,"Rhat"])

system("shutdown -h now")

