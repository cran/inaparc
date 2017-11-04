\name{spaeth}
\alias{spaeth}
\title{
Initialization of cluster prototypes using Spaeth's algorithm
}
\description{
Initializes the cluster prototypes using the centroids that are calculated with Spaeth's algorithm (Spaeth, 1977), which is similar to Forgy's algorithm.
}
\usage{
spaeth(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
In this algorithm, each object in the data set is assigned to one of \var{k} clusters in cyclical fashion. The j-\emph{th} \eqn{(j \epsilon {1,2, \cdots, n}}{(j \epsilon {1,2, \cdots, n}}) object is assigned to the \eqn{(j-1 (mod k) + 1)}\eqn{(j-1 (mod\, k)+1)}-\emph{th} cluster. In contrast to Forgy's method, this method is sensitive to order of data (Celebi et al, 2013).
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{avg} with this function because the cluster prototypes are the averages of sampled objects for each feature.}
  \item{call}{a string containing the matched function call that generates the \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Spaeth, H. (1977). Computational experiences with the exchange method: Applied to four commonly used partitioning cluster analysis criteria, \emph{European J of Operational Rsch.}, 1(1):23-31. doi:\url{https://doi.org/10.1016/S0377-2217(77)81005-9}

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. arXiv:\url{https://arxiv.org/pdf/1209.1960.pdf}
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
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
res <- spaeth(iris[,1:4], k=5)
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