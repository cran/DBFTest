\name{cdf.F}
\alias{cdf.F}
\title{Null cumulative distribution function of the DBF statistic}
\description{Returns the probability of observing a DBF statistic value greater than or equal to an observed value. This function is intended for internal use.}
\usage{
cdf.F(dmat,group.labels,a)
}
\arguments{
  \item{dmat}{Distance matrix of size \code{n} x \code{n}}
  \item{group.labels}{Vector of length \code{n} containing group membership of each sampling unit}
  \item{a}{Observed DBF statistic value for which a p-value is required}
  }
\value{The probability of observing a DBF statistic value greater than or equal to \code{a}.}
\author{Christopher Minas}
\references{Minas, C. and Montana, G. (2012). Distance-based analysis of variance: approximate inference. Submitted to Statistical Analysis and Data Mining. arxiv: http://arxiv.org/abs/1205.2417  }
