% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InvTransform.R
\name{InvTurkey}
\alias{InvTurkey}
\title{Performs a Inverse Turkey power tranformation with predefined Lambda}
\usage{
InvTurkey(x, lambda)
}
\arguments{
\item{x}{A vector which needs to be inverse-tranformed}

\item{lambda}{A fixed lambda value for turkey power transformation}
}
\value{
Returns exponential tranformattion if lambda equals zero and a inverse power transormation otherwise
}
\description{
Takes in a vector and a lambda value to perform Inverse Turkey power transfomation
}
