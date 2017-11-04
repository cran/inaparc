\name{crsamp}
\alias{crsamp}
\title{
Initialization of cluster prototypes using the centers of random samples
}
\description{
Initializes the cluster prototypes matrix using the centers of \var{r} data objects. The options for centers are mean and median of the sampled objects in addition to the objects nearest to the mean of the sampled objects.
}
\usage{
crsamp(x, k, r, ctype)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{r}{an integer for the number of objects to be sampled from the data set. If missing, the default value is 2.}
  \item{ctype}{a string for the type of centroids to be computed. The options are \sQuote{avg} for average, \sQuote{med} for median or \sQuote{obj} for the object nearest to the average. The default is \sQuote{obj}.}
}

\details{
Instead of sampling only one random object as the function \code{\link{rsamp}} does, the function \code{arsamp} randomly samples \var{r} data objects, and then computes the average and median of these sampled objects. The nearest data object to the mean of sampled objects is also found. If \code{ctype} is \option{avg} the mean of the sampled \var{r} objects is assigned as the prototype of first cluster. When \code{ctype} is \option{med} the median of the sampled \var{r} objects is assigned as the prototype of first cluster. If the  \code{ctype} is \option{obj}, the nearest object to the mean of sampled \var{r} objects is assigned as the the prototype of first cluster. The same process is repeated for all of the remaining clusters. The logic behind this novel technique is to avoid to select the outliers in the data set which may occur with random sampling for only one object.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string for the type of used centroid to build the cluster prototypes matrix.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
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
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
# Prototypes are the objects nearest to the mean of
# five randomly sampled objects for each cluster
res <- crsamp(iris[,1:4], k=5, r=5, ctype="obj")
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{sampling for prototype selection}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}