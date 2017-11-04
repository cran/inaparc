\name{get.algorithms}
\alias{get.algorithms}
\title{
Get the names of algorithms in \sQuote{inaparc}
}
\description{
Gets the names of initialization algorithms which are available in the package \sQuote{\pkg{inaparc}}.
}
\usage{
get.algorithms(atype="prototype")
}

\arguments{
  \item{atype}{an string for the type of algorithms. The default value is \option{prototype} for the names of algorithms for initialization of cluster prototypes. Use \option{membership} for the names of algorithms for initialization of hard and fuzzy membership degrees.}
}

\value{a vector containing the names of algorithms.}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\seealso{
 \code{\link{inaparc-package}}
}

\examples{
get.algorithms(atype="prototype")
get.algorithms(atype="membership")
}
