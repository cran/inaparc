\name{ksegments}
\alias{ksegments}
\title{
Initialization of cluster prototypes using the centers of \var{k} segments
}
\description{
Initializes the cluster prototypes matrix using the centers of \var{k} segments (subsets) of the data set. 
}
\usage{
ksegments(x, k, ctype)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{ctype}{a string for the type of centroid. The options are \sQuote{avg} for average and \sQuote{med} for median of the objects in the segments. The default is \sQuote{avg}.}
}

\details{
The first segment consists of the first \code{n/k} objects. The second segment consists of \code{n/k} objects starting from the \code{n/k+1}-\emph{th} object. The process is repeated for \var{k} segments. The centers of \var{k} segments are assigned as the cluster prototypes.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid. Its value is \sQuote{avg} for average or \sQuote{med} for median of the objects in the segments.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}
\seealso{
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{inofrep}},
 \code{\link{inscsf}},
 \code{\link{insdev}},
 \code{\link{kkz}},
 \code{\link{kmpp}},
 \code{\link{ksteps}},
 \code{\link{lastk}},
 \code{\link{lhsrandom}},
 \code{\link{lhsmaximin}},
 \code{\link{mscseek}},
 \code{\link{rsegment}},
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

# Generate the prototypes matrix using the means of segments
res <- ksegments(x=iris[,1:4], k=5, ctype="avg")
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}