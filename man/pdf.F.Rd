\name{pdf.F}
\alias{pdf.F}
\title{Null probability denisty function of the DBF statistic}
\description{Returns the PDF of the DBF statistic under the null hypothesis over a given support.}
\usage{pdf.F(dmat, group.labels,xs)}
\arguments{
  \item{dmat}{Distance matrix of size \code{n} x \code{n}}
  \item{group.labels}{Vector of length \code{n} containing group membership of each sampling unit}
  \item{xs}{Vector of values at which PDF is required}
  }

\details{Analytic expressions for the mean, variance and skewness of the exact permutation distribution 
of the DBF statistic are computed, and the PDF is approximated by the Pearson Type III distribution with these mean, variance and skewness quantities. 
}
\value{A vector of the same length as \code{xs} giving the PDF values at those points.}
\references{Minas, C. and Montana, G. (2012). Distance-based analysis of variance: approximate inference. Submitted to Statistical Analysis and Data Mining. arxiv: http://arxiv.org/abs/1205.2417  }
\author{Christopher Minas}

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

# plot of PDF of DBF: 

xs <- seq(-0.005,0.07,len=1000)	
pdf.dbf <- pdf.F(dmat, group.labels,xs)
plot(xs,pdf.dbf,type="l",ylab="PDF(x)",xlab="x",main="PDF of DBF")

# compare PDF of DBF with sampling distribution using 1000 Monte Carlo permutations:

gmat <- create.G(n,dmat)
H <- scale(make.H(n,group.labels),center=TRUE,scale=FALSE)
dbf.perms <- sapply(1:1000, function(i) {perms <- sample(1:n,replace=FALSE)
							DBF(G=gmat[perms,perms],H)})
hist(dbf.perms,prob=TRUE,add=TRUE)

}
