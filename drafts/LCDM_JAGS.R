library(StanDCM)
library(rstan)
library(R2jags)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

mod.LCDM <- StanLCDM.run(Qmatrix = Qmatrix,
                         response.matrix = respMatrix,
                         iter = 500, warmup = 100, chain.num = 5)

mod.jags <- '

'

##############################################################################
 #              Latent Class Attribute Profile Generator                      #
 ##############################################################################

 generateLatentClasses <- function(qmatrix){
      num_att = length(qmatrix[1,])
      max_att = 2^num_att
      latent_classes = matrix (data=NA, max_att, num_att)
      m <- max_att
      for (a in 1:num_att) {
          m = m/2     # Number of repititions of entries 0 or 1 in one cycle
          anf <- 1
          while (anf < max_att) {
            latent_classes[anf : (anf + m -1),a] <- 0
            anf <- anf + m
            latent_classes[anf : (anf + m -1),a] <- 1
            anf <- anf + m
          }
      }
      rownames(latent_classes) = paste("c", 1:max_att, sep= "")
      latent_classes
 } #### end of function generateLatentClasses

## all.patterns of q entries
gapp <- function(q_matrix){
  K <- ncol(q_matrix)
  q.entries <- as.list(1:K)
  for (k in 1:K) {
    q.entries[[k]] <- sort(unique(c(0, q_matrix[,k])))
  }
  attr.patt <- as.matrix(expand.grid(q.entries))
}
all.patterns <- gapp(Qmatrix)
# q_pattern <- all.patterns
## helper function to create the indicator function of lambda
helper <- function(q_pattern) {
  c(1, vec, )
}

# Useful functions ------------------------------------------------------------
dec2bin = function(decimal_number, nattributes, basevector){
  dec = decimal_number
  profile = matrix(NA, nrow = 1, ncol = nattributes)
  for (i in nattributes:1){
    profile[1,i] =  dec%%basevector[i]
    dec = (dec-dec%%basevector[i])/basevector[i]
  }
  return(profile)
}

bin2dec = function(binary_vector, nattributes, basevector){
  dec = 0
  for (i in nattributes:1){
    dec = dec + binary_vector[i]*(basevector[i]^(nattributes-i))
  }
  return(dec)
}

DINA <- function(){
  for (n in 1:N) {
    for (i in 1:I) {
      for (k in 1:k) {
        w[n,i,k] <- pow(alpha[n,k], Q[i,k]) # weigths
      }
      eta[n,i] <- prod(w[n, i, 1:K])
      p[n,i] <- pow((1 - s[i]), eta[n,i]) * pow(g[i], (1 - eta[n,j]))) # pi LL function
      Y[n,i] <- dbern(p[n,i])
    }
    for (k in 1:k) {
      alpha[n,k] <- all.patterns[c[n], k]
    }
    c[n] ~ dcat(pai[1:C])
  }
  pai[1:C] ~ ddirch(delta[1:C])
  for (i in 1:I) {
    s[i] ~ dbeta(1,1)
    g[i] ~ dbeta(1,1) %_% T( , 1 - s[i])
  }

  ## the posterior predictive model checking
  for (n in 1:N) {
    for (i in 1:I) {
      teststat[n,j] <- pow(Y[n,i] - p[n,j], 2) / (p[n,j] * (1-p[n,j]))
      Y_rep[n,i] ~ dbern(p[n,i])
      teststat_rep[n,i] <- pow(Y_rep[n,j]-p[n,j], 2) / (p[n,j] * (1-p[n,j]))
    }
  }
  teststatsum <- sum(teststat[1:N, 1:I])
  teststatsum_rep <- sum(teststat_rep[1:N, 1:I])
  ppp <- step(teststatsum_rep - teststatsum)
}


LCDM <- function() {
    for (n in 1:N) {
      p[n,1] <- lambda0 + lambda1_1
      p[n,2] <- lambda0 + lambda1_1
      p[n,4] <- lambda0 + lambda1_2
      p[n,5] <- lambda0 + lambda1_2
      p[n,6] <- lambda0 + lambda1_3
      p[n,7] <- lambda0 + lambda1_1 + lambda1_2
      p[n,8] <- lambda0 + lambda1_2 + lambda1_3
      p[n,9] <- lambda0 + lambda1_1 + lambda1_3
      for (i in 1:I) {
        Y[n,i] <- dbern(p[n,i])
        for (k in 1:k) {
          w[n,i,k] <- pow(alpha[n,k], Q[i,k]) # weigths
          pi[i, k] = pi[i, k] + l[eff, index]
        }
      eta[n,i] <- prod(w[n, i, 1:K])
      }
    }


    for (k in 1:k) {
      alpha[n,k] <- all.patterns[c[n], k]
    }
    c[n] ~ dcat(pai[1:C])
}



Qmatrix

p[1,1] = lamba0 + lambda1_1
p[1,2] = lamba0 + lambda1_1
p[1,3] = lamba0 + lambda1_2
p[1,7] = lambda0 + lambda1_2 + lambda1_3 + lambda2_23







