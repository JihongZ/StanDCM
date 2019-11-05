library(StanDCM)

mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 2000, warmup = 1000)
mod.LCDM@stanmodel


## Arguments
stan.model = mod.LCDM
response.matrix = respMatrix
n.sim = 500
n.burnin = 100

StanDCM.ppmc(mod.LCDM, respMatrix, n.sim = 500, n.burnin = 100)
StanDCM.ppmc(mod.LCDM, respMatrix, n.sim = 500, n.burnin = 100, type = "sumscores")
StanDCM.ppmc(mod.LCDM, respMatrix, n.sim = 500, n.burnin = 100, type = "chisq")

posterior.matrix <- as.matrix(stan.model)
Ni <- ncol(response.matrix) # number of item (9)
Np <- nrow(response.matrix) # number of people (1000)
## Expect probability (PI)
PImat.name <- grep(names(stan.model), pattern = "PImat", value = T)
Nc <-  length(PImat.name) / Ni # number of response pattern
Vc.name <- grep(names(stan.model), pattern = "Vc", value = T) # Probability of each response pattern
n.iter = stan.model@stan_args[[1]]$iter # number of iteration
Vc.posterior.matrix <- posterior.matrix[(n.burnin +1):n.iter, Vc.name] # burn out first iteration
PImat.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, PImat.name]

message("Posterior Predictive Model Checking is running: ")

# simulate response matrix based on PImat and Vc
pseudo.sumscore.dist.matrix = NULL
if (is.na(n.sim)) {
  n.sim = Np
}
draw.timer = 0
# create progress bar
# pb <- txtProgressBar(min = 0, max = n.sim, style = 3)
# extract the sumscore frequnency for each draw
pseudo.sumscore.extract <- function() {
  draw.timer <<- draw.timer + 1
  iter <- sample(1:(n.iter-n.burnin), 1, replace = TRUE)
  # select Predicted Probability Matrix
  PImat.select <- matrix(PImat.posterior.matrix[iter,], Ni, Nc)

  # draw from dirichlet distribution and generate class for each person
  Vc.draw <- Vc.posterior.matrix[iter,]
  ## Generate each person's latent class
  pseudo.class.vector<-apply(rmultinom(Np, 1, Vc.draw), 2, function(x){which(x==1)})
  ## each person's probability in each item
  pseudo.response.matrix <- t(PImat.select[, pseudo.class.vector])
  ## create binary response
  pseudo.response.matrix <- t(apply(pseudo.response.matrix, 1, function(x) rbinom(Ni, 1, x)))
  ## Sum all response up for each person
  pseudo.sumscore.vector <- apply(pseudo.response.matrix, 1, sum)
  pseudo.sumscore.dist.vector <- data.frame(table(pseudo.sumscore.vector))
  pseudo.sumscore.dist.vector$time = draw.timer
  pseudo.sumscore.dist.vector
}
# extract bivariate cell frequency
bivariate.chisq.extrat <- function(){
  draw.timer <<- draw.timer + 1
  iter <- sample(1:(n.iter-n.burnin), 1, replace = TRUE)
  # select Predicted Probability Matrix
  PImat.select <- matrix(PImat.posterior.matrix[iter,], Ni, Nc)

  # draw from dirichlet distribution and generate class for each person
  Vc.draw <- Vc.posterior.matrix[iter,]
  ## Generate each person's latent class
  pseudo.class.vector<-apply(rmultinom(Np, 1, Vc.draw), 2, function(x){which(x==1)})
  ## each person's probability in each item
  pseudo.response.matrix <- t(PImat.select[, pseudo.class.vector])
  ## create binary response
  pseudo.response.matrix <- t(apply(pseudo.response.matrix, 1, function(x) rbinom(Ni, 1, x)))

  ## Create 2x2 Contingency Table
  grid_table = t(combn(1:Ni, 2))
  chisq_vec <- sapply(1:nrow(grid_table), function(x){
    i1 <- grid_table[x,1]
    i2 <- grid_table[x,2]
    expected_contingency = c(table(pseudo.response.matrix[,i1],
                                   pseudo.response.matrix[,i2]))
    observed_contingency = c(table(respMatrix[,i1], respMatrix[,i2]))
    chisq = sum((observed_contingency - expected_contingency)^2 / expected_contingency)
    return(chisq)
    }
  )
  bivariate.chisq.dat <- data.frame(cbind(grid_table, chisq_vec))
  colnames(bivariate.chisq.dat) <- c("X1", "X2", "Chisq")
  bivariate.chisq.dat$time = draw.timer
  bivariate.chisq.dat
}

# (2) Set up the style for progree bar
pboptions(type = "txt", style = 3, char = "=")
#（1）extract the pseudo sumscore for n.sim times
pseudo.sumscore.dist.matrix = pbreplicate(n.sim, bivariate.freq.extrat(), simplify = FALSE)
pseudo.sumscore.dist.matrix.long <- do.call(rbind,pseudo.sumscore.dist.matrix)
pseudo.sumscore.dist.matrix.wide <- spread(key = pseudo.sumscore.vector, value = Freq , pseudo.sumscore.dist.matrix.long)
pseudo.sumscore.dist.df <- pseudo.sumscore.dist.matrix.wide[,-1]
#(2) Compute the observed raw score distribution
rep.quantile <- apply(pseudo.sumscore.dist.df, 2, quantile,na.rm = TRUE, probs = c(0.05, 0.5, 0.95))
rep.quantile <- data.frame(t(rep.quantile))
colnames(rep.quantile) <- c("Q5", "Q50", "Q95")
rep.quantile$score <- as.numeric(rownames(rep.quantile))
response.sumscore.obs <- apply(response.matrix, 1, sum)
observed <- data.frame(table(response.sumscore.obs))
colnames(observed) <- c("score", "n.obs")
observed$score <- as.numeric(as.character(observed$score))
# (3) combine them together
ppmc.result <- left_join(observed, rep.quantile, by = "score")
ppmc.result$'WithinPosterior(95%)' <- ifelse(ppmc.result$n.obs > ppmc.result$Q95 | ppmc.result$n.obs < ppmc.result$Q5, FALSE, TRUE)
