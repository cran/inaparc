\name{scseek}
\alias{scseek}
\title{
Initialization of cluster prototypes using SCS algorithm
}
\description{
Initializes the cluster prototypes matrix with the Simple Cluster Seeking (SCS) algorithm (Tou & Gonzales, 1974).
}
\usage{
scseek(x, k, tv)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
  \item{tv}{a number to be used as the threshold distance which is directly input by the user. Also it is possible to compute \var{T}, a threshold distance value with the following options of \code{tv} argument:
  \itemize{
    \item \var{T} is the mean of differences between the consecutive pairs of objects with the option \option{cd1}.  
    \item \var{T} is the minimum of differences between the consecutive pairs of objects with the option \option{cd2}.     
    \item \var{T} is the mean of Euclidean distances between the consecutive pairs of objects divided into \var{k} with the option \option{md}. This is the default if \code{tv} is not supplied by the user.
    \item \var{T} is the range of maximum and minimum of Euclidean distances between the consecutive pairs of objects divided into \var{k} with the option \option{mm}.
  }
 }
}

\details{
The algorithm Simple Cluster Seeking (SCS) (Tou & Gonzales, 1974) is similar to Ball and Hall's algorithm (Ball & Hall, 1967) with an exception for selection of the first object (Celebi et al, 2013). In SCS, the first object in the data set is selected as the prototype of the first cluster. Then, the next object whose distance to the first prototype is greater than \var{T}, a threshold distance value is seeked and assigned as the second cluster prototype, if found. Afterwards, the next object whose distance to already determined prototypes is greater than \var{T} is searched and assigned as the third cluster prototype. The selection process is repeated for determining the prototypes of remaining clusters in similar way. 

Because SCS is sensitive to the order of the data (Celebi et al, 2013), it may not yield good initializations with the sorted data. On the other hand, the distance between the cluster prototypes can be controlled \var{T}, which is an arbitrary number specified by the user. But the problem is that how the user decides on this threshold value. As a solution to this problem in the function \code{scseek}, some internally computed distance measures can be used. (See the section\sQuote{Arguments} above for the available options.)}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix of the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototype matrix contains the objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}
\references{
Ball, G.H. & Hall, D.J. (1967). A clustering technique for summarizing multivariate data, \emph{Systems Res. & Behavioral Sci.}, 12 (2): 153-155.

Tou, J.T. & Gonzalez,R.C. (1974). \emph{Pattern Recognition Principles}. Addison-Wesley, Reading, MA. <ISBN:9780201075861>

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
 \code{\link{scseek2}},
 \code{\link{spaeth}},
 \code{\link{ssamp}},
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
# Run with the threshold value of 0.5
res <- scseek(x=iris[,1:4], k=5, tv=0.5)
v1 <- res$v
print(v1)

# Run with the internally computed default threshold value 
res <- scseek(x=iris[,1:4], k=5)
v2 <- res$v
print(v2)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}