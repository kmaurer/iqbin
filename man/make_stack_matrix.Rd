% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IterativeQuantileSupportFunctions.R
\name{make_stack_matrix}
\alias{make_stack_matrix}
\title{Stack Matrix builder}
\usage{
make_stack_matrix(N, M)
}
\arguments{
\item{N}{number of rows}

\item{M}{number of duplicates}
}
\value{
bin index for new observation
}
\description{
function to make a J matrix for duplicating the N rows of a matrix M times each. (support funciton for iqbin function)
}
