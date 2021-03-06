% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StanDINA.script.R
\name{StanDINA.script}
\alias{StanDINA.script}
\title{A function to generate leave-one-out cross-validation for Stan Model}
\usage{
StanDINA.script(Qmatrix, save.path = getwd(), save.name = "DINA_uninf")
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
The StanLCDM.loofit Function to automate Stan code geneartion for LCDMs with binary resposnes
}
\author{
{Zhehan Jiang, University of Alabama, \email{zjiang17@ua.edu}}
}
