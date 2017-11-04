\name{inofrep}
\alias{inofrep}
\title{
Initialization of cluster prototypes using Inofrep algorithm
}
\description{
Initializes cluster prototypes using Inofrep which is a novel prototypes initialization algorithm using the peaks of frequency polygon of a selected feature.
}
\usage{
inofrep(x, k, sfidx, sfpm, binrule, nbins, tcmethod, tc)
}

\arguments{
  \item{x}{a numeric vector, data frame or matrix.}
  \item{k}{an integer for the number of clusters.}
  \item{sfidx}{an integer specifying the column index of a selected feature which is used for determination of protoypes. If missing, it is internally determined by comparing the peak counts of all features in the data set, and the feature having maximum number of peaks is used as the selected feature.}
  \item{sfpm}{a numeric two-column matrix containing the middle values and frequencies of the peaks of the selected feature, respectively.}
  \item{binrule}{a string containing the name of binning rule to generate the classes of frequency polygons of features in the data set. If missing, \sQuote{sqr} rule is used as the default, and square root of the row number of data matrix is assigned as the number of classes to generate frequency polygons.}
  \item{nbins}{an integer for the number of classes of frequency polygons of features in the data set. It should be given if the binning rule \sQuote{usr} is selected as the threshold computing method. If missing, it is internally assigned by the binning rule given as the input.}
  \item{tcmethod}{a string representing the threshold value computing method which is used to remove small peaks and empty classes. If missing, the defult method is 'min' which assigns the threshold value to the minimum frequency of the classes in a frequency polygon.}
  \item{tc}{a numeric threshold value for removing the small peaks and empty classes. If missing, it is assigned internally by the used threshold computing method if it is described or 1 if it is not described.}
}

\details{
Inofrep, \dfn{initialization on the frequency polygon} of a selected feature is a data dependent semi-deterministic initialization algorithm to improve the computational efficiency in prototype-based hard and fuzzy clustering. In the descriptive statistics, frequency polygons serve the structural information about the data. Since a cluster is a dense region of objects that is surrounded by a region of low density (Tan et al, 2006), the  peaks of a frequency polygon occur in the center of dense regions of data (Aitnouri et al, 1999). Based on this assumption, the algorithm Inofrep uses that the peak values in frequency polygons as the estimates of central tendency locations or the centres of different dense regions, namely the clusters in the data set. Thus, the peak values can be used as the prototypes of clusters.
}

\value{an object of class \sQuote{inaparc}, which is a list consists of the following items:
  \item{v}{a numeric matrix containing the initial cluster prototypes.}
  \item{sfidx}{an integer for the column index of the selected feature, which used for determination of cluster prototypes.}
  \item{ctype}{a string for the type of centroid, which used for assigning the cluster prototypes.}
  \item{call}{a string containing the matched function call that generates this \sQuote{inaparc} object.}
}

\references{
Aitnouri E.M., Wang, S., Ziou, D., Vaillancourt, J. & Gagnon, L. (1999). An algorithm for determination of the number of modes for pdf estimation of multi-modal histograms, in \emph{Proc. of Vision Interface '99}, Trois-Rivieres, Canada, May 1999, p. 368-374. \cr \url{http://www.crim.ca/perso/langis.gagnon/articles/vi99_aitnouri.pdf} 

Tan, P. N., Steinbach, M., & Kumar, V. (2006). Cluster analysis: Basic concepts and algorithms. In \emph{Introduction to Data Mining}, Pearson Addison Wesley. \url{http://www-users.cs.umn.edu/~kumar/dmbook/ch8.pdf}
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\note{
In order to supply the peak matrices directly, the functions \code{\link[kpeaks]{findpolypeaks}} and \code{\link[kpeaks]{rmshoulders}} of the package \sQuote{\pkg{kpeaks}} can be used.
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{forgy}},
 \code{\link{hartiganwong}},
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
 \code{\link{scseek2}},
 \code{\link{ssamp}},
 \code{\link{topbottom}},
 \code{\link{uniquek}},
 \code{\link{ursamp}}
}

\examples{
data(iris)
# set 2nd feature as the selected feature
sfidx <- 2

# generate frequency polygon for the selected feature with user-defined class number
hvals <- kpeaks::genpolygon(iris[,sfidx], binrule="usr", nbins=20)

# Call findpolypeaks for calculating the peaks matrix for the peaks of frequency polygon
resfpp <- kpeaks::findpolypeaks(hvals$mids, hvals$freqs, tcmethod="min")
sfpm <- resfpp$pm

# Call inofrep with the peaks matrix calculated in previous step
v <- inofrep(x=iris[,1:4], k=5, sfidx=sfidx, sfpm=sfpm)$v
print(v)
}

\concept{initialization of cluster prototypes}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}