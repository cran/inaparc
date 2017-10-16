\name{firstk}
\alias{firstk}
\title{
Initialization of cluster prototypes using the first k objects
}
\description{
Initializes the cluster prototypes matrix using the first \var{k} objects at the top of data set. 
}
\usage{
firstk(x, k)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
The technique so-called the \dfn{first method of MacQueen} (MacQueen, 1967) that simply selects the first \var{k} objects as the initial centroids. It is sensitive to the order of data (Celebi et al, 2013). If the data set is already sorted in any order it may result with no good initial prototypes because the data objects are close to each other in a sorted data set. Therefore, shuffling of the data set as a pre-processing step may improve the quality with this initialization technique.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of used centroid to build prototype matrix. Its value is \sQuote{obj} with this function because it returns the selected objects.}
  \item{call}{a string containing the matched function call that generates the object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
MacQueen, J.B. (1967). Some methods for classification and analysis of multivariate observations, in \emph{Proc. of 5-th Berkeley Symp. on Mathematical Statistics and Probability}, Berkeley, University of California Press, 1: 281-297.\cr <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.8619&rep=rep1&type=pdf>

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. <arXiv:https://arxiv.org/pdf/1209.1960.pdf>
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
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
res <- firstk(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}