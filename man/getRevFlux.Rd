\name{getRevFlux}
\alias{getRevFlux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getRevFlux
}
\description{
Given flux of irreversible model the function finds foreward minus backward flux}
\usage{
getRevFlux(model, modirrev, fdirrev)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{   An object of class \code{\link{modelorg}}.}
  \item{modirrev}{ An object of class \code{\link{modelorg}} with only irreversible reactions.
 }
  
  \item{fdirrev}{
fluxes of irreversible model
}
}

\value{return fluxes according to the reversible model
}

\author{Abdelmoneim Amer Desouki}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mod2irrev}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, modirrev, fdirrev) 
{
    fluxes = NULL
    for (r in (react_id(model))) {
        if (!react_rev(model)[react_id(model) == r]) {
            fluxes = rbind(fluxes, cbind(rxn = r, fwd = fdirrev[which(react_id(modirrev) == 
                r)], bwd = 0))
        }
        else {
            fluxes = rbind(fluxes, cbind(rxn = r, fwd = fdirrev[which(react_id(modirrev) == 
                paste(r, "_f", sep = ""))], bwd = fdirrev[which(react_id(modirrev) == 
                paste(r, "_b", sep = ""))]))
        }
    }
    return(fluxes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
