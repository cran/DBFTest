\name{create.G}
\alias{create.G}
\title{Gower's centered inner product matrix}
\description{Returns Gower's centered inner product matrix arising from a given distance matrix.}
\usage{create.G(n, D)}
\arguments{
  \item{n}{Number of sampling units}
  \item{D}{Distance matrix of size \code{n} x \code{n}}
}
\value{Gower's centered inner product matrix of size \code{n} x \code{n}.}

\references{Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika}
\author{Christopher Minas}
