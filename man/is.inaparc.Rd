\name{is.inaparc}
\alias{is.inaparc}
\title{
Checking the object class for \sQuote{inaparc}
}
\description{
Checks whether the given object is an instance of the \code{inaparc} class. 
}
\usage{
is.inaparc(x)
}

\arguments{
  \item{x}{an object to check.}
}

\value{\code{TRUE} if \code{x} is a valid \code{inaparc} object and \code{FALSE} for the other type of object classes.
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{crsamp}},
 \code{\link{inofrep}},
 \code{\link{inscsf}},
 \code{\link{insdev}},
 \code{\link{kkz}},
 \code{\link{kmpp}},
 \code{\link{ksegments}},
 \code{\link{ksteps}},
 \code{\link{lastk}},
 \code{\link{lhsrandom}},
 \code{\link{lhsmaximin}},
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
res <- firstk(x=iris[,1:4], k=5)
is.inaparc(res)

x <- c(1,5,8)
is.inaparc(x)
}
