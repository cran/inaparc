\name{inscsf}
\alias{inscsf}
\title{
Initialization cluster prototypes using Inscsf algorithm
}
\description{
Initializes cluster prototypes with \code{Inscsf} which is a novel prototype initialization algorithm using a selected central tendency measure of a selected feature. For reducing the computational complexity and increasing the accuracy in initialization, the algorithm works on only one feature which can be selected according to its importance in clustering. Furthermore, with a selection mechanism using the distribution of data the algorithm also automatically decides what type of center measure should be used.
}
 
\usage{
inscsf(x, k, sfidx, ctype)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer specifying the number of clusters.}
  \item{sfidx}{an integer specifying the column index of the selected feature. If missing, it is internally determined by comparing the number of unique values for all the features in the data set. The feature having the maximum number of unique values is used as the selected feature.}
  \item{ctype}{a string for the type of the selected center. The options are \sQuote{avg} for average, \sQuote{med} for median or \sQuote{mod} for mode. The default value is \sQuote{avg}.}
}

\details{
The \code{inscsf} is based on a technique so-called "\dfn{initialization using a selected center of a selected feature}". It resembles Ball and Hall's method (Ball and Hall, 1967) for assignment of the first cluster prototype but it differs by the use of two different interval values (\eqn{R_1}{R1} and \eqn{R_2}{R2}) instead of using only one fixed threshold value (\emph{T}) for determining the prototypes of remaining clusters. The technique \code{inscsf} does not require to sort the data set. \eqn{R_1}{R_1} is an interval which is calculated by dividing the distance between the center and maximum of the selected feature (\eqn{x_f}{x_f}) by half of the number of clusters minus 1. 

\deqn{R_1=\frac{max(x_f)-center(x_f)}{(c-1)/2}}{R_1=\frac{max(x_f)-center(x_f)}{(c-1)/2}} 

Similarly, \eqn{R_2}{R_2} is an interval which is calculated by dividing the distance between the maximum and center of the selected feature by half of the number of clusters minus 1. 

\deqn{R_2=\frac{center(x_f)-min(x_f)}{(k-1)/2}}{R_2=\frac{center(x_f)-min(x_f)}{(k-1)/2}}  

These two intervals become equal to each other if the selected feature is normally distributed, and thus, cluster prototypes are located in equidistant positions from each other in the \var{p}-dimensional space of \var{n} data objects.

Depending on the distribution of selected feature, the mean, median or mode of the selected feature can be used to determine the prototype of first cluster. If the type of center measure is not input by the user, it is internally determined according to the distribution of data. Then, the nearest data instance to the center of the selected feature is searched on the selected feature column, and assigned as the prototype of first cluster. 

\deqn{v_1=x_i ;  i = row \;index \;of \;the \;nearest  \;data \;object \;to \;center(x_f))}{v_1=x_i ;  i = row index of the nearest data instance to center(x_f))}

The prototype of an even-numbered cluster is determined by adding the value of \eqn{R_1}{R_1} times the cluster index minus 1 to the first cluster prototype.

\deqn{v_j={(x_{(i+(j-1)}\;R_1)}}{v_i={(x_{(i+(j-1)} R1)}}

On the other hand, \eqn{R_2}{R_2} is used to calculate the prototypes for the odd-numbered clusters. 
\deqn{v_j={(x_{(i+(j-1)}\; R_2)}}{v_i={(x_{(i+(j-1)} R_2)}}
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{sfidx}{an integer for the column index of the selected feature.}
  \item{ctype}{a string for the type of centroid. It is \sQuote{obj} with this function because the prototypes matrix contain contains the selected objects.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\references{
Ball, G.H. & Hall, D.J. (1967). A clustering technique for summarizing multivariate data, \emph{Systems Res. & Behavioral Sci.}, 12 (2): 153-155.

Cebeci, Z., Sahin, M. & Cebeci, C. (2018). Data dependent techniques for initialization of cluster prototypes in partitioning cluster analysis. In Proc. of \emph{4th International Conference on Engineering and Natural Science}, Kiev, Ukraine, May 2018. pp. 12-22.
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\note{
The selected feature can be determined in several ways. The feature with highest number of peaks among the others can be also utilized as the selected feature with this function. For determination of it, the function \code{\link[kpeaks]{findpolypeaks}} of the package \sQuote{\pkg{kpeaks}} can be used.
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{forgy}},
 \code{\link{hartiganwong}},
 \code{\link{inofrep}},
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
 \code{\link{ursamp}}
}

\examples{
data(iris)
# Use the 4th feature as the selected feature
v1 <- inscsf(x=iris[,1:4], k=5, sfidx=4)$v
print(v1)

# Use the internally selected feature
v2 <- inscsf(x=iris[,1:4], k=5)$v
print(v2)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}