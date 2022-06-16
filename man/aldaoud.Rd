\name{aldaoud}
\alias{aldaoud}
\title{
Initialization of cluster prototypes using Al-Daoud's algorithm
}
\description{
Initializes the cluster prototypes matrix using the variance-based algorithm proposed by Al-Daoud (Al-Daoud, 2005). 
}
\usage{
aldaoud(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
At first, the algorithm finds the feature having the greatest variance and sorts the data set on this feature in any order. Then it divides the data set into \code{n/k}-length \var{k} segments. The medians of the segments are assigned as the protoypes of clusters. Al-Daoud's algorithm is likely to be effective only for data sets in which the variability is mostly on one dimension because it considers only one feature with the highest variance (Celebi et al, 2013).
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{sfidx}{an integer for the column index of the feature with the highest variance.}
  \item{ctype}{a string for the type of used centroid to determine the cluster prototypes. It is \sQuote{med} with this function.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Al-Daoud, M.B. (2005). A new algorithm for cluster initialization, in \emph{Proc. of 2nd World Enformatika Conf.}, pp.74-76.

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. arXiv:\url{https://arxiv.org/pdf/1209.1960.pdf}
}

\seealso{
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
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
res <- aldaoud(iris[,1:4], k=5)
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