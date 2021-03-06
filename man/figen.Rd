\name{figen}
\alias{figen}
\title{
Initialization of membership degrees over class range of a selected feature
}
\description{
Initializes the membership degrees matrix by using the class range of the coefficient of variation of a selected feature in the data set being processed.
}
\usage{
figen(x, k, mtype, sfidx)
}

\arguments{
  \item{x}{an data.frame or matrix for the data set.}
  \item{k}{an integer for the number of clusters.}
  \item{mtype}{a character representing the type of membership degrees to be generated. The default type is \option{f} for generating fuzzy membership matrix. Use \option{h} for creating an hard (crisp) membership matrix.}
  \item{sfidx}{an integer for the column index of a selected feature. The default is the column index of a feature whose coefficient of variation is the maximum among all features in the data set.}
}

\details{
The function \code{figen} generates a numeric matrix containing the fuzzy initial membership degrees.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{u}{a numeric matrix containing the initial membership degrees.}
  \item{sfidx}{an integer for the column index of the selected feature, which used for random sampling.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Cebeci, Z. (2018), "Initialization of Membership Degree Matrix for Fast Convergence of Fuzzy C-Means Clustering", In Proc. of \emph{2018 International Conference on Artificial Intelligence and Data Processing (IDAP)}, IEEE, Sep. 2018, pp. 1-5., doi: \href{https:/doi.org/10.1109/IDAP.2018.8620920}{10.1109/IDAP.2018.8620920}
}

\seealso{
 \code{\link{imembrand}},
 \code{\link{imembones}}
}

\examples{
data(iris)

# Generate a fuzzy membership matrix using the 1st feature
u <- figen(iris[,1:4], k=5, sfidx=1)$u
head(u)
tail(u)

# Generate a fuzzy membership matrix using the internally determined feature
res <- figen(iris[,1:4], k=5)
u <- res$u
head(u)
tail(u)
}

\concept{initialization of membership degrees}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}