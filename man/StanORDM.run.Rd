% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StanORDM.run.R
\name{StanORDM.run}
\alias{StanORDM.run}
\title{Generate Stan code and Run the estimation for ORDM}
\usage{
StanORDM.run(
  Qmatrix,
  response.matrix,
  script.path = NA,
  save.path = getwd(),
  save.name = "ORDM_uninf",
  iter = 1000,
  warmup = 0,
  chain.num = 3,
  init.list = "random",
  control.list = NA
)
}
\arguments{
\item{Qmatrix}{the Q-matrix specified for the LCDM}

\item{response.matrix}{the response matrix}

\item{script.path}{the path to save the .stan file to somewhere; the default path is getwd()}

\item{save.path}{the path to save the .stan file to somewhere; the default path is getwd()}

\item{save.name}{name of the .stan}

\item{iter}{number of iterations}

\item{warmup}{number of warmup iterations}

\item{chain.num}{number of MCMC chains. Default is 3.}

\item{init.list}{initial variables. Default is random. Other options include cdm}

\item{control.list}{constrains put on the multigroup ORDM model.}
}
\value{
a. stan file saved at the specified path
}
\description{
The StanORDM.run Function to automate Stan code geneartion for LCDMs with binary resposnes
}
\author{
{Zhehan Jiang, University of Alabama, \email{zjiang17@ua.edu}}
}
