\name{getGpr1iso}
\alias{getGpr1iso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose the the smallest isoenzyme
}
\description{
preprocessing of model to get the minimal cost iso enzyme and eliminate the OR from GPR's.
}
\usage{
getGpr1iso(model, cpx_stoich, MW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  object of class \code{\link{modelorg}} representing the metabolic network.
}
  \item{cpx_stoich}{
 giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 1/10/2015
               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
              'stoich': number of subunits of each gene in the same order in 'genes'
}
  \item{MW}{
MW measurement for gene using readfaa.r, \code{\link{calc_MW}} or gene annotation [units g/mmol]
  }
}
\value{
return a LIST:
  \item{gpr1iso }{data frame containing the choosen term (complex or isoenzyme)}
  \item{mod_cpx_mw }{molecular weight of isoenzyme}
  \item{rgst}{matrix of reaction gene with stoichiometry}
}
\author{
Abdelmoneim Amer Desouki
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{getccFBA_mat}}
}
