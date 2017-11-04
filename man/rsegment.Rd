\name{rsegment}
\alias{rsegment}
\title{
Initialization of cluster prototypes using a randomly selected segment
}
\description{
Initializes the cluster prototypes matrix using using a \var{k}-length segment of data set consists of consecutive objects that starts with a randomly sampled data object.
}
\usage{
rsegment(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The function \code{rsegment} randomly samples one data object as the prototype of first cluster, and then it assigns the next \var{k-1} data objects as the prototype of remaining clusters.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototype matrix contains the sampled objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{forgy}},
 \code{\link{hartiganwong}},
 \code{\link{inofrep}},
 \code{\link{inscsf}},
 \code{\link{insdev}},
 \code{\link{kkz}},
 \code{\link{kmpp}},
 \code{\link{ksegments}},
 \code{\link{ksteps}},
 \code{\link{lastk}},
 \code{\link{lhsmaximin}},
 \code{\link{lhsrandom}},
 \code{\link{maximin}},
 \code{\link{mscseek}},
 \code{\link{rsamp}},
 \code{\link{scseek}},
 \code{\link{scseek2}},
 \code{\link{spaeth}},
 \code{\link{ssamp}},
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
res <- rsegment(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}