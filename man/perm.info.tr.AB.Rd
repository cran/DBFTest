\name{perm.info.tr.AB}
\alias{perm.info.tr.AB}
\title{Function to compute the mean, variance and skewness of the exact permutation distribution the statistic \code{T=trace(AB)}}
\description{Returns the permutational mean, variance and skewness of the statistic \code{T=trace(AB)}, where \code{A} and \code{B} are square, symmetric and centered (both column- and row-centered) matrices of the same size. This is achieved by using analytic formulae that arise by considering the exact permutation distribution of \code{T}. The exact permutation 
distribution refers to the discrete distribution obtained by recomputing \code{T} for all possible simultaneous permutations of the rows and columns of \code{B}. The function also returns
\code{trace(B)}.This function is intended for internal use. 
}

\usage{
perm.info.tr.AB(A, B, n)
}
\arguments{
  \item{A}{Square, symmetric and centered matrix of size \code{n} x \code{n}}
  \item{B}{Square, symmetric and centered matrix of size \code{n} x \code{n}}
  \item{n}{Number of rows and columns in \code{A} and \code{B}}
}
\value{
	\item{mean.T}{Permutational mean of \code{T}}
	\item{variance.T}{Permutational variance of \code{T}}
	\item{skewness.T}{Permutational skewness of \code{T}}
	\item{tr.G}{trace(\code{B})}	
}

\references{Kazi-Aoual, F., Hitier, S., Sabatier, R. and Lebreton, J. D. (1995) Refined approximations to permutation tests for multivariate inference. Computational Statistics and Data Analysis. Elsevier Science B.V.}
\author{Christopher Minas}

