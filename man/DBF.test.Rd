\name{DBF.test}
\alias{DBF.test}
\title{DBF Test}
\description{Returns the results of applying the DBF test of no difference between groups to an observed distance matrix with group membership information.}
\usage{DBF.test(dmat, group.labels, no.permutations)}
\arguments{
  \item{dmat}{Distance matrix of size \code{n} x \code{n}}
  \item{group.labels}{Vector of length \code{n} containing group membership of each sampling unit}
  \item{no.permutations}{Number of permutations to use in order to obtain a p-value. If 0, an approximation to the exact permutation distribution is used.}
}
\details{If the number of permutations is set to 0, the exact permutation distribution of the DBF statistic which would be obtained by permuting the rows and columns of \code{dmat}, and 
recomputing the DBF statistic for all \code{n}! permutations, is considered in order to obtain a p-value. Using analytic expressions for the permutational mean, variance and skewness of this exact 
permutation distribution, the exact permutation distribution is approximated by the continuous Pearson type III distribution. The observed statistic is compared to this continuous 
distribution, and the approximate p-value obtained. If permutations are used instead, the discrete Monte Carlo permutation distribution is obtained with the given number of permutations.
As before, the observed statistic is compared to this, and the permutation p-value obtained.}
\value{A vector with elements \code{dbf.statistic} and \code{dbf.p.value}.}
\references{Minas, C., Waddell, S. J. and Montana, G. (2011). Distance-based differential analysis of gene curves. Bioinformatics

Minas, C. and Montana, G. (2012). Distance-based analysis of variance: approximate inference. Submitted to Statistical Analysis and Data Mining. arxiv: http://arxiv.org/abs/1205.2417}
\author{Christopher Minas}
\note{The distance matrix used must be symmetric and contain non-negative values.}
\examples{
# set seed for reproducibility:

set.seed(101)

# create multivariate normal dataset of size 99 x 10:

n <- 99
p <- 10
X <- matrix(rnorm(n*p, mean=0, sd=1),nrow=n)	

# set 3 equal-sized groups:

n1 <- n2 <- n3 <- n/3
group.labels <- c(rep(1,n1),rep(2,n2),rep(3,n3))	

# define Manhattan distance matrix:

dmat <- as.matrix(dist(X,"manhattan"))

# apply DBF test of equality between groups using approximation:

DBF.test(dmat,group.labels,no.permutations=0)
#dbf.statistic   dbf.p.value 
#   0.02198219    0.40518124 

# apply DBF test of equality between groups using 1000 permutations:

DBF.test(dmat,group.labels,no.permutations=1000)
#dbf.statistic   dbf.p.value 
#   0.02198219    0.41700000 
}

