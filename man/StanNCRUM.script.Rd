% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StanNCRUM.script.R
\name{StanNCRUM.script}
\alias{StanNCRUM.script}
\title{Generate Stan code and Run the estimation for ORDM}
\usage{
StanNCRUM.script(Qmatrix, save.path = getwd(), save.name = "NCRUM_uninf")
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
