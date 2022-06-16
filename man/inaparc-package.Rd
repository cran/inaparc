\name{inaparc-package}
\alias{inaparc-package}
\docType{package}
\title{Initialization Algorithms for Partitioning Cluster Analysis}
\description{Partitioning clustering algorithms divide data sets into k subsets or partitions which are so-called clusters. They require some initialization procedures for starting to partition the data sets. Initialization of cluster prototypes is one of such kind of procedures for most of the partitioning algorithms. Cluster prototypes are the data elements, i.e. centroids or medoids, representing the clusters in a data set. In order to initialize cluster prototypes, the package \sQuote{\pkg{inaparc}} contains a set of the functions that are the implementations of widely-used algorithms in addition to some novel techniques. Initialization of fuzzy membership degrees matrices is another important task for starting the probabilistic and possibilistic partitioning algorithms. In order to initialize membership degrees matrices required by these algorithms, the package \sQuote{\pkg{inaparc}} contains a number of functions for most of the data independent and dependent initialization techniques (Borgelt, 2005) which are categorized as the linear time-complexity and loglinear time complexity-initialization methods in Celebi et al (2013). 
}

\details{Clustering is one of the most widely used exploratory statistical analysis in data mining. Its goal is to explore the groups of objects that are similar to each other within the group but different from the objects in other groups. According to a common taxonomy, the existing clustering algorithms are classified in two groups: Hierarchical and Non-hierarchical (or flat) algorithms (Rokah & Maimon, 2005). As a dominant subfamily of non-hierarchical algorithms, the partitioning clustering algorithms divide data objects into a pre-defined number of clusters, which are the non-overlapping subsets of data. Although the choice of an appropriate algorithm for any clustering task depends on many criteria or purposes. When data size and dimensions are the concerned criteria, the non-hierarchical algorithms may be more practical way of clustering the large size and high dimensional data sets because they quickly process the large data sets when compared to the hierarchical clustering algorithms.

As the most crowded group of the partitioning clustering tools, the prototype-based algorithms partition data objects into clusters in which each data object is more similar to its prototype than the prototypes of other clusters. On clustering context, a prototype is a typical data item that represents or characterizes a cluster (Tan et al. 2006). Usually, it can be regarded as the most central data point in a data subspace so-called cluster. The prototype of a cluster is so often a centroid, i.e., the mean of all the objects in a cluster. On the other hand, centroids can not be computed for non-numeric data, i.e., on nominal or ordinal data. In such case, medoids can be used as the prototypes of clusters (Tan et al, 2006). 

Initialization or seeding is a process for selecting the starting values of cluster prototypes matrix which serves the initial representatives of clusters. It is an important task in partitioning cluster analysis because it is known that the final clustering result is to be highly sensitive to the initial prototypes of the clusters (Khan, 2012). When the prototypes are chosen to be equal or close to the actual centers of clusters in a data set, the partitioning converges quickly and yields quality results. Contrarily, poor initializations of prototype matrix may result with no-good quality of final partitions. 

In fuzzy and possibilistic clustering, an object is a member of all clusters in varying degrees of membership instead of being a member of only one cluster. A membership degrees matrix is required by the fuzzy clustering algorithms, i.e., Fuzzy C-means (FCM) (Bezdek, 1981). Initialization of membership degrees for starting FCM and its various variants must satisfy the following constraints:  

\deqn{u_{ij}\in[0,1]; 1\le i \le n, 1\le j \le k}{u_{ij}\in[0,1]; 1\le i \le n, 1\le j \le k}
\deqn{\sum\limits_{j=1}^k u_{ij}=1; 1\le i \le n}{\sum\limits_{j=1}^k u_{ij}=1; 1\le i \le n}
\deqn{0<\sum\limits_{i=1}^n u_{ij} < n ; 1\le j \le k}{0<\sum\limits_{i=1}^n u_{ij} < n ; 1\le j \le k}

Membership degrees matrices are usually initialized with the techniques based on random number generating as the function \code{\link{imembrand}} does. In addition to these common techiques, a novel technique using the information from synthetically produced classes over a selected feature is provided in the package \sQuote{\pkg{inaparc}}. The novel technique which is implemented in \code{\link{figen}} may contribute to the fast convergence of the clustering algorithms when compared to the random sampling based techniques. The package also serves the functions for building hard or crisp membership degrees which can be used for testing purposes.
}

\author{Zeynel Cebeci, Cagatay Cebeci}

\references{
Bezdek J.C. (1981). Pattern recognition with fuzzy objective function algorithms. Plenum, NY, 256 p. <ISBN:0306406713>

Borgelt, C., (2005). \emph{Prototype-based classification and clustering}. Habilitationsschrift zur Erlangung der Venia legendi fuer Informatik, vorgelegt der Fakultaet fuer Informatik der Otto-von-Guericke-Universitaet Magdeburg, Magdeburg, 22 June 2005. url:\url{https://borgelt.net/habil/pbcc.pdf}

Cebeci, Z. (2018). Initialization of Membership Degree Matrix for Fast Convergence of Fuzzy C-Means Clustering", In Proc. of \emph{2018 International Conference on Artificial Intelligence and Data Processing (IDAP)}, IEEE, Sep. 2018, pp. 1-5., \doi{https://doi.org/10.1109/IDAP.2018.8620920}

Cebeci, Z., Sahin, M. & Cebeci, C. (2018). Data dependent techniques for initialization of cluster prototypes in partitioning cluster analysis. In Proc. of \emph{4th International Conference on Engineering and Natural Science}, Kiev, Ukraine, May 2018. pp. 12-22. 

Rokah, L. & Maimon, O. (2005). Clustering methods. In \emph{Data Mining and Knowledge Discovery Handbook (ed. O. Maimon)}, Springer US. pp. 321-352. \doi{https://doi.org/10.1.1.149.9326}

Tan, P. N., Steinbach, M., & Kumar, V. (2006). Cluster analysis: Basic concepts and algorithms. In \emph{Introduction to Data Mining}. Pearson Addison Wesley. url:\url{https://www-users.cse.umn.edu/~kumar/dmbook/ch8.pdf}

Khan, F. (2012). An initial seed selection algorithm for k-means clustering of georeferenced data to improve replicability of cluster assignments for mapping application. \emph{Applied Soft Computing}, 12 (11) : 3698-3700. \doi{https://doi.org/10.1016/j.asoc.2012.07.021}

Celebi, M.E., Kingravi, H.A. & Vela, P.A. (2013). A comparative study of efficient initialization methods for the K-means clustering algorithm, \emph{Expert Systems with Applications}, 40 (1): 200-210. arXiv:\url{https://arxiv.org/pdf/1209.1960.pdf}
}

\seealso{
 \code{\link{aldaoud}},
 \code{\link{ballhall}},
 \code{\link{crsamp}},
 \code{\link{firstk}},
 \code{\link{forgy}},
 \code{\link{hartiganwong}},
 \code{\link{imembones}},
 \code{\link{imembrand}},
 \code{\link{figen}},
 \code{\link{inofrep}},
 \code{\link{inscsf}},
 \code{\link{insdev}},
 \code{\link{is.inaparc}},
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

\concept{initialization of cluster prototypes}
\concept{cluster seeding techniques}
\concept{initialization of membership degrees matrix}
\concept{prototype-based clustering}
\concept{partitioning clustering}
\concept{partitional clustering}
\concept{non-hierarchial clustering}
\concept{flast clustering}
\concept{cluster analysis}
\concept{unsupervised learning}

\keyword{cluster}
