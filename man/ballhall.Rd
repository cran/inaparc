\name{ballhall}
\alias{ballhall}
\title{
Initialization of cluster prototypes using Ball & Hall's algorithm
}
\description{
Initializes the prototypes of clusters by using the cluster seeding algorithm which has been proposed by Ball & Hall (1967).
}
\usage{
ballhall(x, k, tv)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{tv}{a number to be used as \var{T}, a threshold distance value. It is directly input by the user. Also it is possible to compute \var{T} with the following options of \code{tv} argument:
  \itemize{
    \item \var{T} is the mean of differences between the consecutive pairs of objects with the option \option{cd1}.  
    \item \var{T} is the minimum of differences between the consecutive pairs of objects with the option \option{cd2}.     \item \var{T} is the mean of Euclidean distances between the consecutive pairs of objects divided into \var{k} with the option \option{md}. This is the default if \code{tv} is not supplied by the user.
    \item \var{T} is the range of maximum and minimum of Euclidean distances between the consecutive pairs of objects divided into \var{k} with the option \option{mm}.
  }
 }
}

\details{
In the Ball and Hall's algorithm (Ball & Hall, 1967), the center of gravity of data is assigned as the prototype of first cluster. It then passes the data objects in arbitrary order and takes an object as the next prototype if it is \var{T} units far from the previously selected prototypes. The purpose of using \var{T}, the distance threshold, is to make the cluster protoypes at least \var{T} units away from each other. Ball & Hall's method may be sensitive to the order of data, and moreover, deciding for an appropriate value of \var{T} is is also difficult (Celebi et al, 2013). As the solutions to this problem, the function \code{ballhall} in this package computes a \var{T} value using some distance measures, if it is not specified by the user (for details, see the section \sQuote{Arguments} above.)
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string for the type of used centroid. It is \sQuote{obj} with this function because the created cluster prototypes matrix contains the selected objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\references{
Ball, G.H. & Hall, D.J. (1967). A clustering technique for summarizing multivariate data, \emph{Systems Res. & Behavioral Sci.}, 12 (2): 153-155.

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. <arXiv:https://arxiv.org/pdf/1209.1960.pdf>
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{aldaoud}},
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
 \code{\link{ursamp}},
}

\examples{
data(iris)
# Run with a user described threshold value
v1 <- ballhall(x=iris[,1:4], k=5, tv=0.6)$v
print(v1)

# Run with the internally computed default threshold value
v2 <- ballhall(x=iris[,1:4], k=5)$v
print(v2)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}