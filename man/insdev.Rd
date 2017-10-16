\name{insdev}
\alias{insdev}
\title{
Initialization of cluster prototypes using Insdev algorithm
}
\description{
\code{Insdev} is a novel algorithm that initializes the cluster prototypes by using the standard deviation of a selected feature. The selected feature is the most important feature in regard of variation. For this purpose the coefficients of variation of the features are compared, and then the feature with highest coefficient of variation is selected for further processes.
}
\usage{
insdev(x, k, sfidx)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{sfidx}{an integer specifying the column index of the selected feature. Here, in this function we use the feature with high variability as the selected feature because it dominates the clustering results (Khan, 2912). If missing, so it is internally determined by comparing the coefficents of variation for all the features in the data set. The feature having the maximum coefficient of variation is used as the selected feature.}
}

\details{
At first the algorithm computes the mean of the selected feature (\eqn{\bar{x_{s}}}{\bar{x_{s}}}) and then seeks the object whose distance is minimum to \eqn{\bar{x_{s}}}{\bar{x_{s}}} as the prototype of first cluster. The prototypes of remaining clusters are determined by using a stepping range (\emph{R}), computed from the standard deviation of selected feature with the formula \eqn{R=1/2\sigma_{x_{s}}/k}{R=1/2\sigma_{x_{s}}/k}. The prototype of second cluster is the object whose distance is minimum to \eqn{\bar{x_{s}} + (i-1)\; R}{\bar{x_{s}} + (i-1) R}, where \emph{i} is the cluster index. The prototype of third cluster is the object whose distance is minimum to \eqn{\bar{x_{s}} - i \; R}{\bar{x_{s}} - i R} in the opposite direction to previous prototype. The prototypes remaining clusters are cyclically determined in similar way.

Since it produces the same prototypes in each run of it, \code{insdev} is a deterministic algorithm. Therefore, this characteristic of the algorithm provides replicability in  initialization procedure.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{sfidx}{an integer for the column index of the selected feature, used in the calculations.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because the cluster prototypes are the objects sampled from the data set.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Khan, F. (2012). An initial seed selection algorithm for k-means clustering of georeferenced data to improve replicability of cluster assignments for mapping application. \emph{Applied Soft Computing}, 12 (11) : 3698-3700. <doi:https://doi.org/10.1016/j.asoc.2012.07.021>
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
res <- insdev(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}