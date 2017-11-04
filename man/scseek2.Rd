\name{scseek2}
\alias{scseek2}
\title{
Initialization of cluster prototypes using SCS algorithm over a selected feature
}
\description{
Initializes the cluster prototypes matrix with the Simple Cluster Seeking (SCS) algorithm (Tou & Gonzales, 1974) over a selected feature.
}
\usage{
scseek2(x, k, sfidx, tv)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
  \item{sfidx}{an integer specifying the column index of the selected feature for random sampling. If missing, it is internally determined by comparing the coefficients of variation of all features in the data set. The feature having the maximum coefficent of variation is used as the selected feature.}
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
The \code{scseek2} is a novel variant of the function \code{\link{scseek}} based on the Simple Cluster Seeking (SCS) algorithm (Tou & Gonzales, 1974). It differs from SCS that the distances and threshold value are computed over a selected feature having the maximum coefficient of variation, instead of using all the features.  
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix of the initial cluster prototypes.}
  \item{sfidx}{an integer for the column index of the selected feature, which used for random sampling.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototype matrix contains the sampled objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{proclus} object.}
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
 \code{\link{mscseek}},
 \code{\link{rsamp}},
 \code{\link{rsegment}},
 \code{\link{scseek}},
 \code{\link{spaeth}},
 \code{\link{ssamp}},
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
# Run over 4th feature with the threshold value of 0.5
res <- scseek2(x=iris[,1:4], k=5, sfidx=4, tv=0.5)
v1 <- res$v
print(v1)

# Run with the internally computed default threshold value 
res <- scseek2(x=iris[,1:4], k=5)
v2 <- res$v
print(v2)

}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}