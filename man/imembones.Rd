\name{imembones}
\alias{imembones}
\title{
Initialization of a crisp membership matrix using a selected cluster
}
\description{
Initializes a crisp membership degrees matrix which is used to start a partitional clustering algorithm.
}
\usage{
imembones(n, k, mtype, numseed)
}

\arguments{
  \item{n}{an integer for the number of objects in the data set.}
  \item{k}{an integer for the number of clusters.}
  \item{mtype}{a string representing the type of crisp initialization for a selected cluster. The default is 'hrc'. The alternatives are 'hfc' in which all objects are assumed as the member of the first cluster, and 'hlc' in which all objects are assumed as the member of the last cluster.}
  \item{numseed}{a number to be used for the seed of RNG.}
}

\details{
The function \code{imembones} generates a numeric membership degrees matrix containing the crisp initial values for a selected cluster.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{u}{a numeric matrix containing the crisp initial membership degrees of the objects to \var{k} clusters.}
  \item{sfidx}{an integer for the column index of the selected feature, which used for random sampling.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{imembrand}},
 \code{\link{figen}}
}

\examples{
# Generate membership degrees matrix whose last column contains crisp
# membership degrees
u <- imembones(n=10, k=5, mtype="hlc")$u
head(u)
tail(u)

# Generate membership degrees matrix using a seed number
u <- imembones(n=10, k=5, mtype="hrc", numseed=123)$u
head(u)
tail(u)
}

\concept{initialization of membership degrees}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}