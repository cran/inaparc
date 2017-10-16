\name{forgy}
\alias{forgy}
\title{
Initialization of cluster prototypes using Forgy's algorithm
}
\description{
Initializes the cluster prototypes using the centers that are calculated with Forgy's algorithm (Forgy, 1965), which is the earliest algorithm for seeding the clusters in the standard K-means clustering. 
}
\usage{
forgy(x, k)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
In this algorithm, each object in the data set is randomly assigned to one of \var{k} clusters, and then the mean of the objects assigned to the clusters are used as the initial cluster prototypes. The algorithm lacks of theoretical basis, and the clusters generated randomly may have no internal homogeneity (Celebi et al, 2013).
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{avg} with this function because the cluster prototypes are the averages of sampled objects for each cluster.}
  \item{call}{a string containing the matched function call that generates the object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Forgy, E.W. (1965). Cluster analysis of multivariate data: Efficiency vs interpretability of classification, \emph{Biometrics}, 21 (3) : 768-769.

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. <arXiv:https://arxiv.org/pdf/1209.1960.pdf>
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
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
res <- forgy(iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{sampling for prototype selection}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}