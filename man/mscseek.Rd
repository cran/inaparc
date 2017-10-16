\name{mscseek}
\alias{mscseek}
\title{
Initialization of cluster prototypes using the modified SCS algorithm
}
\description{
Initializes the cluster prototypes matrix using a modified version of the Simple Cluster Seeking (SCS) algorithm proposed by Tou & Gonzales(1974). While SCS uses a fixed threshold distance value \var{T} for selecting all candidates of clusters, modified SCS recomputes \var{T} with the average Euclidean distances between the previously determined prototypes. This adjustment makes possible to select more cluster prototypes when compared to SCS.
}
\usage{
mscseek(x, k, tv)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
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
This is a modification of the Simple Cluster Seeking (SCS) algorithm (Tou & Gonzalez, 1974). The algorithm selects the first object in the data set as the prototype of the first cluster. Then, next object whose distance to the first prototype is greater than a threshold distance value is searched and assigned as the second cluster prototype. Instead of using a fixed the \var{T}, threshold distance value as SCS does, the modified SCS recomputes the \var{T} by the average Euclidean distances between the previously determined prototypes of clusters. The next object whose distance to the previously selected object is greater than the adjusted \var{T} is searched and assigned as the third cluster prototype. The selection process is repeated for the remaining clusters in similar way. The method is sensitive to the order of the data, it may not yield good initializations with the ordered data. 
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
Tou, J.T. & Gonzalez,R.C. (1974). \emph{Pattern Recognition Principles}. Addison-Wesley, Reading, MA. <ISBN:9780201075861>
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
# Run with the threshold value of 0.1
res <- mscseek(x=iris[,1:4], k=5, tv=0.1)
v1 <- res$v
print(v1)

# Run with the internally computed default threshold value 
res <- mscseek(x=iris[,1:4], k=5)
v2 <- res$v
print(v2)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}