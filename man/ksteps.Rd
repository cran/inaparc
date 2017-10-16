\name{ksteps}
\alias{ksteps}
\title{
Initialization of cluster prototypes using the centers of \var{k} blocks
}
\description{
Initializes the cluster prototypes matrix using the centers of objects in \var{k} blocks that are generated with a kind of systematic sampling method as described in the section \sQuote{Details}. 
}
\usage{
ksteps(x, k, ctype)
}

\arguments{
  \item{x}{a numeric data frame or matrix to be processed.}
  \item{k}{an integer for the number of clusters.}
  \item{ctype}{a string for the type of centroid. The options are \sQuote{avg} for average and \sQuote{med} for median of the objects in the blocks. The default is \sQuote{avg}.}
}

\details{
The algorithm \code{ksteps} is similar to \code{\link{ksegments}} but it differs for the selection of the members of the segments or blocks. The objects whose row indexes are \code{1, 1+k, 1+2k,...} are assigned to the first segment, and then the objects whose row indexes are \code{2, 2+k, 2+2k,...} to the second block. In this way, \var{k} blocks of the objects are formed. The centers of these \var{k} blocks are assigned as the cluster prototypes.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix.}
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
res <- ksteps(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}