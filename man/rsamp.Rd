\name{rsamp}
\alias{rsamp}
\title{
Initialization of cluster prototypes using simple random sampling
}
\description{
Initializes the cluster prototypes matrix using the randomly selected \var{k} objects from the data set.
}
\usage{
rsamp(x, k)
}

\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
}

\details{
The function \code{rsamp} generates a protoype matrix using the \var{k} objects which are randomly sampled from the data set without replacement. Simple random sampling (SRS), also so-called the \dfn{second method of MacQueen} in clustering context, assumes that cluster areas have a high density; in consequence, the good candidates of the cluster prototypes can be sampled from these dense regions of data with a higher chance (Celebi et al, 2013). SRS is probably the most common approach to initialize prototype matrices. So, it can be seen a \emph{de facto standard} because it has been widely applied with the basic K-means algorithm for the years. Since SRS has no rule to avoid to select the outliers or the objects close to each other, it may result with no good initializations. Before initialization of SRS, multivariate outliers removal on the data set as a data pre-processing step may be helpful to avoid for selection of the outliers, but increases the computational cost.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{ctype}{a string representing the type of centroid, which used to build prototype matrix. Its value is \sQuote{obj} with this function because it samples the objects only.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
MacQueen, J.B. (1967). Some methods for classification and analysis of multivariate observations, in \emph{Proc. of 5-th Berkeley Symp. on Mathematical Statistics and Probability}, Berkeley, University of California Press, 1: 281-297.\cr <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.308.8619&rep=rep1&type=pdf>

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. <arXiv:https://arxiv.org/pdf/1209.1960.pdf>
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
res <- rsamp(x=iris[,1:4], k=5)
v <- res$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitional clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{Statistics|cluster}