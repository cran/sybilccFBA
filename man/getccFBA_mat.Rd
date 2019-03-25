\name{getccFBA_mat}
\alias{getccFBA_mat}

\title{
getccFBA_mat: get ccFBA model as a constraint matrix
}
\description{
given  Kcat vector (in 1/Sec), and molecular weights and optionally protein complex
stoichiometry formulate a ccFBA problem(Maximize biomass). Improvements over MOMENT: 
1- use multifunctioning enzyme constraint. 2-use one enzyme (the cheapest) for each 
reaction catalayzed my multiple isoenzymes and thus eliminate under estimation of 
rules of form ((A AND B) OR (A AND C)) 3-optionally include protein complex stoichiometry
}
\usage{
getccFBA_mat(model, mod2 = NULL, Kcat, MW = NULL, selected_rxns = NULL, verboseMode = 2,
objVal = NULL, RHS = NULL, solver = SYBIL_SETTINGS("SOLVER"), medval = NULL, 
cpx_stoich = NULL, C_mu_coef = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{   An object of class \code{\link{modelorg}}.}
  \item{mod2}{ An object of class \code{\link{modelorg}} with only irreversible reactions.
        It can be sent  to save time of recalculating it with each call.}
  \item{Kcat}{ kcat values in unit 1/S. Contains three slots: reaction id,direction(dirxn),value(val)}
  \item{MW}{ list of molecular weights of all genes, using function calc_MW, in units g/mol}

  \item{selected_rxns}{optional parameter used to select a set of reactions not all, list of react_id}
  \item{verboseMode}{
    An integer value indicating the amount of output to stdout:
    0: nothing, 1: status messages, 2: like 1 plus with more details,
    3: generates files of the LP problem.\cr
    Default: \code{2}.
    }
  \item{RHS}{ the budget C, for EColi 0.27}
 \item{objVal}{when not null the problem will be to find the minimum budget that give the specified 
      objective value(biomass)}

 \item{solver}{
       Single character string giving the solver package to use.  See
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("SOLVER")}.
     }
 \item{medval}{ median of Kcat values , used for missing values}

  \item{cpx_stoich}{
       giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 
               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
               'stoich': number of subunits of each gene in the same order in 'genes'
   Default: NULL (e.g. all stoichiometry coefficients equal one)
   }
  \item{C_mu_coef}{
  used to have C as a linear function of mu (biomass) : C = RHS + C_mu_coef*Biomass
  default 0 (e.g. ignored).
  }
}
\details{
variable names:
g1_x_1: indicates that gene x is in one-to-one relation
gr_x_r:indicates that  gene x catalyzes more than one reaction and
 this variable is the portion that catalyzes reaction r.
gmr_x_cnt: variable for a gene catalyzing more than one reaction (i.e cnt reactions)
           and this variable is the sum of individual parts.
   cpx_r: used to represent a complex catalyzing reaction r.        
}
\value{
return a LIST, 
\item{LHS }{constraint matrix, consists of S matrix and ccFBA constraints}
\item{rlb,rub }{bounds for constraints (i.e rows)}
\item{clb,cub }{bounds for variables (i.e columns)}
\item{rxnkcatRow }{mapping of Kcat values to constraints}
\item{rgst }{complex stoichiometry as returned by \code{\link{getGpr1iso}}}
\item{rxnGeneCol }{rxn_id,gene,Col : used to get MCrowding}
\item{geneCol }{low level structure, used for debugging}
\item{mod_objc_ind }{row index of constraint on FBA model objective}
\item{cc_ind }{capacity constraint index}

}
\references{
Desouki, Abdelmoneim. "Algorithms for improving the predictive power of flux balance analysis." PhD diss., 2016.
}
\author{
Abdelmoneim Amer Desouki
}
\seealso{
 \code{\link{getGpr1iso}}, \code{\link{cfba_moment_mr}},\code{\link{simulate_EColi}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ccFBA }% use one of  RShowDoc("KEYWORDS")
