\name{hartiganwong}
\alias{hartiganwong}
\title{
Initialization of cluster prototypes using Hartigan-Wong's algorithm
}
\description{
Initializes the cluster prototypes matrix using the Hartigan-Wong's algorithm (Hartigan & Wong, 1979). 
}
\usage{
hartiganwong(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
Firstly, the algorithm computes the center of gravity of data and the distances of data objects to this center. Then, it sorts the data set in any order of the computed distances. The prototypes of \var{k} clusters are determined by using the formula (\eqn{1 + (i-1) n/k)}{1 + (i-1) (n/k)}), where \var{i} and \var{n} stand for the index of a cluster and the number of data rows, respectively. This algorithm leads to increase in the computational cost due to complexity of sorting, which is \eqn{O(n\;log(n))}{O(n log(n))} (Celebi et al, 2013).
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string for the type of used centroid to determine the cluster prototypes. It is \sQuote{obj} with this function because the generated prototype matrix contains the selected objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Hartigan, J.A. & Wong, W.A., (1979). Algorithm AS 136: A K-means clustering algorithm, \emph{J of the Royal Statistical Society}, C 28 (1): 100-108.  

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. arXiv:\url{https://arxiv.org/pdf/1209.1960.pdf}
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{forgy}},
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
res <- hartiganwong(iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{sampling for prototype selection}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}