\name{imembrand}
\alias{imembrand}
\title{
Initialization of membership matrix using simple random sampling
}
\description{
Initializes the membership degrees matrix which is used to start a fuzzy and possibilistic partitioning clustering algorithm.
}
\usage{
imembrand(n, k, mtype, numseed)
}

\arguments{
  \item{n}{an integer for the number of objects in the data set.}
  \item{k}{an integer for the number of clusters.}
  \item{mtype}{a string for any of three random initialization methods. The default method is \option{f1} for fuzzy memberships. The options are \option{f2} and \option{f3} for fuzzy memberships and \option{h} for hard (crisp) memberships.}
  \item{numseed}{a number to be used for the seed of RNG.}
}

\details{
The function \code{imembrand} generates a numeric matrix containing the initial membership degrees by using simple random sampling technique.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{u}{a numeric matrix containing the crisp initial membership degrees of \var{n} objects to \var{k} clusters.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{imembucr}},
 \code{\link{imembones}}
}

\examples{
data(iris)
n <- dim(iris)[1]

# Generate a fuzzy membership degrees matrix using default values
u <- imembrand(n=n, k=5)$u
head(u)
tail(u)

# Generate a fuzzy membership degrees matrix using the method 3
u <- imembrand(n=n, k=5, mtype="f3", numseed=123)$u
head(u)
tail(u)

# Generate a crisp membership degrees matrix 
u <- imembrand(n=n, k=5, mtype="h")$u
head(u)
tail(u)
}

\concept{initialization of membership degrees}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}