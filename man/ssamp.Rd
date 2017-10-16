\name{ssamp}
\alias{ssamp}
\title{
Initialization of cluster prototypes using systematic random sampling
}
\description{
Initializes the cluster prototypes matrix using the systemically sampled data objects. 
}
\usage{
ssamp(x, k)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The function \code{ssamp} generates a prototype matrix using the sytematic random sampling technique. Since the data objects are enough away from each other with this technique it may provide better initializations than the simple random sampling. The first object is randomly sampled from the top \code{n/k} objects of data set and assigned as the prototype of first cluster. The prototypes of remaining clusters are the objects whose row indexes are \eqn{v_{1} + i\;(n/k)}{v_{1} + i (n/k)}, where \var{v1} and \var{i} are the index of first selected object and index of cluster, respectively.
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
 \code{\link{rsegment}},
 \code{\link{scseek}},
 \code{\link{scseek2}},
 \code{\link{spaeth}},
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
res <- ssamp(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}