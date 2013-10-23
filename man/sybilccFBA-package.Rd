\name{sybilccFBA-package}
\alias{sybilccFBA-package}
\alias{sybilccFBA}
\docType{package}
\title{
  Cost Constrained Flux Balance Analysis(ccFBA)
}
\description{
The package \code{sybilccFBA} implements some methods to get cost constrained fluxes.
It is required to supply the molecular weights. It can be calculated from genome data using function \code{\link{calc_MW}}.
Also requires kinetic data along with the model.
}
\details{
\tabular{ll}{
Package: \tab sybilccFBA\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2013-06-03\cr
License: \tab GPL Version 3\cr
LazyLoad: \tab yes\cr
Depends: \tab \code{\link{sybil}}, methods\cr
}
}
\author{Abdelmoneim Amer Desouki}
\keyword{ package }
\seealso{
	\code{\link{sybil}} 
	\code{\link{cfba_moment}} 
}