% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COMET_functions.R
\name{fit.CTMC}
\alias{fit.CTMC}
\title{This function optimally fits 3 CTMC models to data (1st phase, 2nd phase,
MET range)}
\usage{
fit.CTMC(data.input, MET.range, opt.cutoff)
}
\arguments{
\item{data.input}{input data to use}

\item{MET.range}{range where MET takes place}

\item{opt.cutoff}{optimal cutoff of highly variable genes}
}
\value{
final trajectories, lambda_E, mu, and lambda_M respectively
}
\description{
This function optimally fits 3 CTMC models to data (1st phase, 2nd phase,
MET range)
}
