\name{lhsmaximin}
\alias{lhsmaximin}
\title{
Initialization of cluster prototypes using Maximin LHS
}
\description{
Initializes the cluster prototypes matrix using the Maximin version of Latin Hypercube Sampling (LHS). A square grid containing possible sample points is a Latin Square (LS) if there is only one sample in each row and each column. LHS is a generalized version of LS, which has been developed to generate a distribution of collections of parameter values from a multidimensional distribution. LHS generates more efficient estimates of desired parameters than simple Monte Carlo sampling (Carnell, 2016). 
}
\usage{
lhsmaximin(x, k, ncp)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{ncp}{an integer determining the number of candidate points used in the search by maximin LHS algorithm.}
}

\details{
LHS aims at initial cluster centers whose coordinates are well spread out in the individual dimensions (Borgelt, 2005). It is the generalization of Latin Square for an arbitrary number of dimensions (features).  When sampling a function of \var{p} features, the range of each feature is divided into \var{k} equally probable intervals. \var{k} samples are then drawn such that a Latin Hypercube is created. 

The current version of the function \code{lhsmaximin} in this package uses the results from the \code{\link[lhs]{maximinLHS}} function from the \sQuote{\pkg{lhs}} library created by Carnell (2016). Once the uniform samples are created by the \code{\link[lhs]{maximinLHS}}, they are transformed to normal distribution samples by using the quantile functions. But all the features in the data set may not be normally distributed, instead they may fit to different distributions. In such cases, the transformation for any feature should be specisific to its distribution. Determination of the distribution types of features is planned in the future versions of the function \sQuote{\code{lhsmaximin}}. 
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string for the type of used centroid to determine the cluster prototypes. It is \sQuote{obj} with this function.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Borgelt, C., (2005). \emph{Prototype-based classification and clustering}. Habilitationsschrift zur Erlangung der Venia legendi fuer Informatik, vorgelegt der Fakultaet fuer Informatik der Otto-von-Guericke-Universitaet Magdeburg, Magdeburg, 22 June 2005. url:\url{http://www.borgelt.net/habil/pbcc.pdf}

Carnell, R., (2016). lhs: Latin Hypercube Samples. R package version 0.14. \url{https://CRAN.R-project.org/package=lhs}
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
res <- lhsmaximin(iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{latin hypercube sampling}
\concept{initialization of cluster prototypes}
\concept{sampling for prototype selection}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}