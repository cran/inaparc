\name{kkz}
\alias{kkz}
\title{
Initialization of cluster prototypes using KKZ algorithm
}
\description{
Initializes the cluster prototypes matrix using \sQuote{KKZ} algorithm proposed by Katsavounidis et al (1994).
}
\usage{
kkz(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
}

\details{
The function \code{kkz} is an implementation of the cluster seeding algorithm which has been proposed by Katsavounidis et al (1994). As the first cluster prototype, the algorithm so-called \sQuote{KKZ} selects one data object on the edges of data. It is the object having the greatest squared Euclidean norm in the function \code{kkz}. The second cluster prototype is the farthest object from the previously selected object. After assignment of the prototypes of first two clusters, the distances of all of the remaining objects to them are computed. The object which is the farthest from its nearest prototype is assigned as the third prototype. The above process is repeated for selecting the prototypes of remaining clusters in the same way. The algorithm \sQuote{KKZ} is considered to be sensitive to the outliers in the data set.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototypes are the selected objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Katsavounidis, I., Kuo, C. & Zhang, Z. (1994). A new initialization technique for generalized Lloyd iteration. \emph{IEEE Signal Processing Letters}, 1 (10): 144-146. url:\url{https://pdfs.semanticscholar.org/0103/d3599757c77f6f3cbe3daf2470f13419cd90.pdf}
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
res <- kkz(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}