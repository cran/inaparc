\name{maximin}
\alias{maximin}
\title{
Initialization of cluster prototypes using Maximin algorithm
}
\description{
Initializes the cluster prototypes matrix by using the Maximin algorithm.
}
\usage{
maximin(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The main idea of the \dfn{Maximin} algorithm is to isolate the cluster prototypes that are farthest apart (Philpot, 2001). The algorithm randomly samples one data object from the data set and assigns it as the first cluster prototype. The prototype of second cluster is determined as the data object which is farthest from the first prototype. Then, the remaining part of data set is scanned for the data objects whose distances are minimum to the previously selected prototypes. The object having the maximum of minimum distances is assigned the prototype of third cluster. The same procedure is repeated for determining the prototypes of other clusters (Spaeth, 1997; Gonzales, 1985; Duda et al, 2000, Celebi et al, 2013).

The algorithm generally works well with circular shaped clusters whose radius are smaller than the separation between clusters.  However, it is very sensitive to the order of object in data sets. Also it is computationally expensive because each time once a new cluster prototype is selected, the distances must be computed for every object from every cluster prototype (Philpot, 2001). In order to contribute to the solutions of this problem, the current implementation of \code{maximin} includes a simple control that if an object has the minimum distance of zero, the seeking procedure is no more continued to compute the distances for the remaining objects. This control may speed the algorithm up with the \code{maximin} function in this package.  
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix of the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototype matrix contains the objects.}
  \item{call}{a string containing the matched function call that generates the object \sQuote{inaparc}.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}
\references{
Spaeth, H. (1977). Computational experiences with the exchange method: Applied to four commonly used partitioning cluster analysis criteria, \emph{European Journal of Operational Research} 1 (1): 23-31. doi:\url{https://doi.org/10.1016/S0377-2217(77)81005-9}

Gonzalez, T. (1985), Clustering to minimize the maximum intercluster distance, \emph{Theoretical Computer Science} 38 (2-3): 293-306. url:\url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.366.8183&rep=rep1&type=pdf}

Duda, R.O., Hart, P.E. & Stork, D.G. (2000). \emph{Pattern Classification}, Wiley-Interscience. <ISBN:978-0-471-05669-0>

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. arXiv:\url{https://arxiv.org/pdf/1209.1960.pdf}

Philpot, W. (2001). Topic 8: Clustering/Unsupervised Classification in \emph{Lecture Notes, CEE 615: Digital Image Processing - Jan 2001, Cornell Univ.}, url:\url{http://ceeserver.cee.cornell.edu/wdp2/cee6150/Monograph/615_08_UnSuper_Class_rev01.pdf}
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
res <-maximin(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}