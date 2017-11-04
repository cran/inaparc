\name{ursamp}
\alias{ursamp}
\title{
Initialization of cluster prototypes using random sampling on each future
}
\description{
Initializes the cluster prototypes matrix by using random uniform sampling for each of \var{p} features in the data set, independently.
}
\usage{
ursamp(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The \code{ursamp} generates the prototypes by binding randomly sampled values for each of \var{p} features, independently. In this novel approach proposed by the authors of the package, an object is randomly sampled from data set and the value of first feature is assigned as the value of first feature of the first prototype. Then next object is sampled and the value of second feature of the sampled object is assigned as the value of second feature of the first prototype. The sampling process is repeated for the other features in similar way. Afterwards the same sampling procedure is repeated for determining the prototypes of remaining clusters.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroids in the prototype matrix. Its value is \sQuote{obj} with this function because it returns objects.}
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
 \code{\link{rsegment}},
 \code{\link{scseek}},
 \code{\link{scseek2}},
 \code{\link{ssamp}},
 \code{\link{spaeth}},
 \code{\link{topbottom}},
 \code{\link{uniquek}}
}

\examples{
data(iris)
res <- ursamp(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}