\name{lhsrandom}
\alias{lhsrandom}
\title{
Initialization of cluster prototypes using random LHS
}
\description{
Initializes the cluster prototypes matrix using the random version of Latin Hypercube Sampling (LHS). A square grid containing possible sample points is a Latin Square (LS) if there is only one sample in each row and each column. LHS is a generalized version of LS, which has been developed to generate a distribution of collections of parameter values from a multidimensional distribution. LHS generates more efficient estimates of desired parameters than simple Monte Carlo sampling (Carnell, 2016). 
}
\usage{
lhsrandom(x, k)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for specifying the number of clusters.}
}

\details{
LHS aims at initial cluster centers whose coordinates are well spread out in the individual dimensions (Borgelt, 2005). LHS is the generalization of Latin Square for an arbitrary number of dimensions (features). When sampling a function of \var{p} features, the range of each feature is divided into \var{k} equally probable intervals. \var{k} samples are then drawn such that a Latin Hypercube is created. 

The current version of the function \code{lhsrandom} in this package uses the results from the \code{\link[lhs]{randomLHS}} function from the  R package \sQuote{\pkg{lhs}} (Carnell, 2016), which contains several variants of LHS. Once the uniform samples are created by the \code{\link[lhs]{randomLHS}}, they are transformed to normal distributed samples by using the quantile functions with \code{\link[stats]{qnorm}}. But all the features in the data set may not be normally distributed, instead they may have the different type of distributions. In such cases, the transformation of any feature should be specific to its distribution. Determination of the distribution types of features is planned in the future versions of the function \sQuote{\code{lhsrandom}}. 
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
 \code{\link{lhsmaximin}},
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
res <- lhsrandom(iris[,1:4], k=5)
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