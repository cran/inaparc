\name{topbottom}
\alias{topbottom}
\title{
Initialization of cluster prototypes using the top and bottom objects
}
\description{
Initializes the cluster prototypes matrix using the alternately selected \var{k} objects from the top and bottom of the data set. 
}
\usage{
topbottom(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The function combines the \code{\link{firstk}} and \code{\link{lastk}} techniques. It takes the first object of the data set as the prototype of first cluster, and then the last object as the prototype of second cluster. This rotating assigment process continues until the prototypes of \var{k} clusters are assigned.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototype matrix contains the objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\note{
If the sorted data set is used, the function \code{topbottom} may yield better initializations when compared to the functions \code{\link{firstk}} and \code{\link{lastk}}.
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
 \code{\link{spaeth}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
res <- topbottom(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}