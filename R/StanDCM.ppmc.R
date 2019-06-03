#' @title A function to calculate pecentage of extreme p-values for sumscore distribution
#'
#' @description
#' The StanDCM.ppmc Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param stan.model A rStan object
#' @param response.matrix the response matrix used by RStan Object
#' @param n.sim number of simulations for Posterior Predictive Model Checking
#' @param n.burnin number of burn-ins
#' @param plot.option logical. whether to provide a plot for ppmc
#' @return p-values tables
#'
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#'
#' @export
#' @examples
#' load("data.RData")
#' Qmatrix<-cbind(Qmatrix,rep(1,9)); Qmatrix[1,1]<-0
#' dim(respMatrix)
#' misspecifiedQmatrix <- Qmatrix
#' misspecifiedQmatrix[1:6,] <- 1-Qmatrix[1:6,]
#' misspecifiedQmatrix[1,3] = 0
#' mod2 <- StanDINA.run(misspecifiedQmatrix,response.matrix = respMatrix,iter=100,init.list='cdm', chain.num = 3, warmup = 20)
#' StanDCM.ppmc(stan.model = mod2, response.matrix = respMatrix, n.sim = 1000, n.burnin = 1, plot.option = FALSE)
#' end - start


StanDCM.ppmc <- function(stan.model, response.matrix, n.sim = 6000, n.burnin = 20, plot.option=FALSE) {
  Install.package(c("Rlab", "MCMCpack", "tidyr", "dplyr", "pbapply"))
  if(plot.option == TRUE) Install.package("ggplot2")

  mod1 <- stan.model
  posterior.matrix <- as.matrix(mod1)
  Ni <- ncol(response.matrix)
  Np <- nrow(response.matrix)
  PImat.name <- grep(names(mod1), pattern = "PImat", value = T)
  Nc <-  length(PImat.name) / Ni
  Vc.name <- grep(names(mod1), pattern = "Vc", value = T)
  n.sim = n.sim
  n.iter = mod1@stan_args[[1]]$iter

  if (n.burnin >= n.iter) {
    stop('Number of iteration less than number of burn-ins')
  }

  Vc.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, Vc.name]
  PImat.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, PImat.name]
  # simulate response matrix based on PImat and Vc
  pseudo.sumscore.dist.matrix = NULL
  if (is.na(n.sim)) {
    n.sim = Np
  }
  time = 0
  # create progress bar
  # pb <- txtProgressBar(min = 0, max = n.sim, style = 3)
  pseudo.sumscore.extract <- function() {
    time <<- time + 1
    iter <- sample(1:(n.iter-n.burnin), 1, replace = TRUE)
    # select Predicted Probability Matrix
    PImat.select <- matrix(PImat.posterior.matrix[iter,], Ni, Nc)

    # draw from dirichlet distribution and generate class for each person
    Vc.draw <- rdirichlet(1, rep(1, Nc))
    pseudo.class.vector<-apply(rmultinom(Np, 1, Vc.draw), 2, function(x){which(x==1)})
    pseudo.response.matrix <- t(PImat.select[, pseudo.class.vector])
    pseudo.response.matrix <- t(apply(pseudo.response.matrix, 1, function(x) rbinom(Ni, 1, x)))
    pseudo.sumscore.vector <- apply(pseudo.response.matrix, 1, sum)
    pseudo.sumscore.dist.vector <- data.frame(table(pseudo.sumscore.vector))
    pseudo.sumscore.dist.vector$time = time
    pseudo.sumscore.dist.vector
  }
  # (2) Set up the style for progree bar
  pboptions(type = "txt", style = 3, char = "=")
  #（1）extract the pseudo sumscore for n.sim times
  pseudo.sumscore.dist.matrix = pbreplicate(n.sim, pseudo.sumscore.extract(), simplify = FALSE)
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
  colnames(observed) <- c("score", "obs")
  observed$score <- as.numeric(as.character(observed$score))
  # (3) combine them together
  df <- left_join(observed, rep.quantile, by = "score")
  df$'WithinPosterior(95%)' <- ifelse(df$obs > df$Q95 | df$obs < df$Q5, FALSE, TRUE)


  if(plot.option == TRUE){
    p <- ggplot(df, aes(x = score + 1, y = obs)) +
          theme_bw() +
          labs(x = "Raw score", y = "Number of examinees", title = "The observed and replicated raw score distributions")

    # (2) Add the entire replicated raw score distributions using violin plot
    # and jittered data points
    p <- p + geom_violin(data = pseudo.sumscore.dist.matrix.long, aes(x = pseudo.sumscore.vector, y = Freq))
    p <- p + geom_jitter(data = pseudo.sumscore.dist.matrix.long, aes(x = pseudo.sumscore.vector, y = Freq, color = pseudo.sumscore.vector),
                         position = position_jitter(width = 0.12), alpha = 0.15) + theme(legend.position = "none")

    # (3) Overlay 5%, 50% (median), and 95% quantiles of each replicated raw
    # score distribution
    p <- p +
      geom_point(aes(y = Q50), size = 4, shape = 2) +
      geom_line(aes(y = Q50),size = 1, linetype = "dashed")
    p <- p + geom_line(aes(y = Q5), size = 1, linetype = "dotted")
    p <- p + geom_line(aes(y = Q95), size = 1, linetype = "dotted")

    # (4) Overlay the observed raw score distribution
    p <- p + geom_point(size = 4) + geom_line(size = 1)
    print(p)
  }

  df
}







