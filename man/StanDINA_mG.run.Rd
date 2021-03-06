% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StanDINA_mG.run.R
\name{StanDINA_mG.run}
\alias{StanDINA_mG.run}
\title{Generate Stan code and Run the estimation for ORDM}
\usage{
StanDINA_mG.run(
  Qmatrix,
  response.matrix,
  GroupID,
  fixeditem.vector = NA,
  class.equal = T,
  script.path = NA,
  save.path = getwd(),
  save.name = "DINA_uninf_multiG",
  iter = 1000,
  warmup = 0,
  chain.num = 3,
  init.list = "random",
  control.list = NA
)
}
\arguments{
\item{Qmatrix}{the Q-matrix specified for the LCDM}

\item{save.path}{save the .stan file to somewhere; the default path is getwd()}

\item{save.name}{name the .stan}
}
\value{
a. stan file saved at the specified path
}
\description{
The StanLCDM.script Function to automate Stan code geneartion for LCDMs with binary resposnes
}
\author{
{Zhehan Jiang, University of Alabama, \email{zjiang17@ua.edu}}
}
