\name{getRevFlux}
\alias{getRevFlux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 get fluxes of reversible model
}
\description{
get fluxes of reversible model from fluxes of irreversible model
}
\usage{
getRevFlux(model, irrxns, fdirrev)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
metabolic network of type \code{\link{modelorg}}
}
  \item{irrxns}{
list of ids of reactions in the irreversible model
}
  \item{fdirrev}{
flux distribution of irreversible reaction
}
}
\value{reaction id in the model and the forward and backward flux, 
  the net flux equals fwd-bwd.
}
\author{
Abdelmoneim Amer Desouki}
