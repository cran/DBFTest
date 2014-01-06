\name{make.H}
\alias{make.H}
\title{Square matrix encoding group membership information}
\description{Returns a matrix encoding group membership information. This function is intended for internal use.}
\arguments{
  \item{n}{Number of sampling units}
  \item{group.labels}{Vector of group labels, of length \code{n}}
}
\value{The square matrix of size \code{n} encoding group membership.}
\author{Christopher Minas}
\references{
Minas, C., Waddell, S. J. and Montana, G. (2011). Distance-based differential analysis of gene curves. Bioinformatics

Minas, C. and Montana, G. (2012). Distance-based analysis of variance: approximate inference. Submitted to Statistical Analysis and Data Mining. arxiv: http://arxiv.org/abs/1205.2417  
}