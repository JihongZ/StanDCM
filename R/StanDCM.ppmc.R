#' @title A function to calculate pecentage of extreme p-values for sumscore distribution
#'
#' @description
#' The StanDCM.ppmc Function to automate Stan code geneartion for LCDMs with binary resposnes
#'
#' @param stan.model A rStan object
#' @param response.matrix the response matrix used by RStan Object
#' @param n.sim number of simulations for Posterior Predictive Model Checking
#' @param n.burnin number of burn-ins
#' @param type The test statistics to perform PPMC. The default is "sumscores". Setting "chisq" will calculate the bivariate item Chi square.
#' @param plot.option logical. whether to provide a plot for ppmc using ggplot2
#' @return p-values tables
#'
#' @author {Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
#' @import Rlab MCMCpack tidyr dplyr pbapply
#' @export
#' @examples
#' \dontrun{
#' load("data.RData")
#' Qmatrix<-cbind(Qmatrix,rep(1,9)); Qmatrix[1,1]<-0
#' dim(respMatrix)
#' misspecifiedQmatrix <- Qmatrix
#' misspecifiedQmatrix[1:6,] <- 1-Qmatrix[1:6,]
#' misspecifiedQmatrix[1,3] = 0
#' mod2 <- StanDINA.run(misspecifiedQmatrix,response.matrix = respMatrix,iter=100,init.list='cdm', chain.num = 3, warmup = 20)
#' StanDCM.ppmc(stan.model = mod2, response.matrix = respMatrix, n.sim = 1000, n.burnin = 1, plot.option = FALSE)
#'}
#'


StanDCM.ppmc <- function(stan.model, response.matrix,
                         n.sim = NULL, n.burnin = NULL,
                         plot.option=FALSE, type = "sumscores") {
  if(plot.option == TRUE) Install.package("ggplot2")

  if (is.null(n.sim)){
    n.sim = stan.model@stan_args[[1]]$iter
  }

  if (is.null(n.burnin)){
    n.burnin = as.integer(n.sim / 10)
  }

  if (n.burnin >= n.sim) {
    stop('Number of simulations less than number of burn-ins')
  }

  posterior.matrix <- as.matrix(stan.model)
  Ni <- ncol(response.matrix)
  Np <- nrow(response.matrix)
  PImat.name <- grep(names(stan.model), pattern = "PImat", value = T)
  Nc <-  length(PImat.name) / Ni
  Vc.name <- grep(names(stan.model), pattern = "Vc", value = T)
  n.iter = stan.model@stan_args[[1]]$iter
  Vc.posterior.matrix <- posterior.matrix[(n.burnin +1): n.iter, Vc.name]
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

  if (type == "sumscores") {
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
    extract_Fun <- pseudo.sumscore.extract

    # (2) Set up the style for progree bar
    pboptions(type = "txt", style = 3, char = "=")
    #（1）extract the pseudo sumscore for n.sim times
    pseudo.sumscore.dist.matrix = pbreplicate(n.sim, extract_Fun(), simplify = FALSE)
    pseudo.sumscore.dist.matrix.long <- do.call(rbind,pseudo.sumscore.dist.matrix)
    pseudo.sumscore.dist.matrix.wide <- spread(key = pseudo.sumscore.vector, value = Freq , pseudo.sumscore.dist.matrix.long, fill = 0)
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

  }else if (type == "chisq"){
    # extract the chisq for each itempair for each draw
    bivariate.chisq.extract <- function(){
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
    extract_Fun <- bivariate.chisq.extract

    # (2) Set up the style for progree bar
    pboptions(type = "txt", style = 3, char = "=")
    #（1）extract the pseudo sumscore for n.sim times
    pseudo.sumscore.dist.matrix = pbreplicate(n.sim, extract_Fun(), simplify = FALSE)
    pseudo.sumscore.dist.matrix.long <- do.call(rbind,pseudo.sumscore.dist.matrix)
    ppmc.result = pseudo.sumscore.dist.matrix.long %>%
      group_by(X1, X2) %>%
      summarise(
        Chisq_post_mean = mean(Chisq),
        Chisq_post_Q5 = quantile(Chisq, prob = 0.05, na.rm = T),
        Chisq_post_Q50 = quantile(Chisq, prob = 0.5, na.rm = T),
        Chisq_post_Q95 = quantile(Chisq, prob = 0.95, na.rm = T),
      )
  }

  if(plot.option == TRUE){
    p <- ggplot(ppmc.result, aes(x = score + 1, y = n.obs)) +
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

  message("Posterior Predictive Model Checking has finished!")
  return(ppmc.result)
}







