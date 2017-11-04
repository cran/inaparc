\name{kmpp}
\alias{kmpp}
\title{
Initialization of cluster prototypes using K-means++ algorithm
}
\description{
Initializes the cluster prototypes matrix by using K-means++ algorithm which has been proposed by Arthur and Vassilvitskii (2007).
}
\usage{
kmpp(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
K-means++ (Arthur & Vassilvitskii, 2007) is usually reported as an efficient approximation algorithm in overcoming the poor clustering problem with the standard K-means algorithm. K-means++ is an algorithm that merges MacQueen's second method with the \sQuote{Maximin} method to initialize the cluster prototypes (Ji et al, 2015). K-means++ initializes the cluster centroids by finding the data objects that are farther away from each other in a probabilistic manner. In K-means++, the first cluster protoype (center) is randomly assigned. The prototypes of remaining clusters are determined with a probability of \eqn{{md(x')}^2/\sum_{k=1}^{n} md({x_k})^2}{{md(x')}^2/\sum_{i=1}^{n} md({x_i})^2}, where \eqn{md(x)}{md(x)} is the minimum distance between a data object and the previously computed prototypes.

The function \code{kmpp} is an implementation of the initialization algorithm of K-means++ that is based on the code\sQuote{k-meansp2.R}, authored by M. Sugiyama. It needs less execution time due to its vectorized distance computations.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototypes are the objects selected by the algorithm.}
  \item{call}{a string containing the matched function call that generates this sQuote{inaparc} object.}
}

\references{
Arthur, D. & Vassilvitskii. S. (2007). K-means++: The advantages of careful seeding, in \emph{Proc. of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms}, p. 1027-1035. url:\url{http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf}

M. Sugiyama, \sQuote{mahito-sugiyama/k-meansp2.R}. url:\url{https://gist.github.com/mahito-sugiyama/ef54a3b17fff4629f106}
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
res <- kmpp(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{k-means++ algorithm}
\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}