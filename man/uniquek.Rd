\name{uniquek}
\alias{uniquek}
\title{
Initialization of cluster prototypes over the unique values
}
\description{
Initializes the cluster prototypes matrix using the randomly sampled data objects over the unique values of a selected feature. 
}
\usage{
uniquek(x, k, sfidx)
}

\arguments{
 \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{sfidx}{an integer specifying the column index of the selected feature for random sampling. If missing, it is internally determined by comparing the number of unique values for all the features in the data set. The feature having the maximum number of unique values is used as the selected feature.}
}

\details{
The set of unique values of the selected feature is determined, and then \var{k} objects were randomly sampled from this set.
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
 \code{\link{ssamp}},
 \code{\link{topbottom}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
# Run with the internally selected feature
res <- uniquek(x=iris[,1:4], k=5)
v <- res$v
print(v)

# Run with the 1st feature
res <- uniquek(x=iris[,1:4], k=5, sfidx=1)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}