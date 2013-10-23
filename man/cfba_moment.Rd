\name{cfba_moment}

\alias{cfba_moment}

\encoding{utf8}

\title{ Function: cfba_moment: implement MOMENT method}
\description{
 This function uses GPR, kcat, and molecular weights to calculate fluxes 
 according to MOMENT method. 
}

\usage{
cfba_moment(model,mod2=NULL, Kcat,MW=NULL,
selected_rxns=NULL,verboseMode=2,objVal=NULL,
RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL)
}

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
obective value(biomass)}

 \item{solver}{
       Single character string giving the solver package to use.  See
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("SOLVER")}.
     }
 \item{medval}{ median of Kcat values , used for missing values}
 
}
\details{
Main steps
    1- Add variables for all genes
    2- for each selected reaction: parse gpr,
    3- Add variables accordingly and constraints
	4-Add solvant constraint

}
\value{
returns a list containing slots:
	 prob:problem object that contains data and model
	 sol: solution of the problem.
     geneCol: mapping of genes to variables in the problem.
}

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{Adadi, R., Volkmer, B., Milo, R., Heinemann, M., & Shlomi, T. (2012). 
Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network 
with Kinetic Parameters, 8(7). doi:10.1371/journal.pcbi.1002575}

\seealso{
 \code{\link{modelorg}},  
\code{\link{optimizeProb}}
}

\examples{
\dontrun{
	library(sybilccFBA)
	data(iAF1260)
	model= iAF1260
 	data(mw)
 	data(kcat)
 	 mod2=mod2irrev(model)
  
	uppbnd(mod2)[react_id(mod2)=="R_EX_glc_e__b"]=1000
	uppbnd(mod2)[react_id(mod2)=="R_EX_glyc_e__b"]=0
	uppbnd(mod2)[react_id(mod2)=="R_EX_ac_e__b"]=0
	uppbnd(mod2)[react_id(mod2)=="R_EX_o2_e__b"]=1000
	lowbnd(mod2)[react_id(mod2)=="R_ATPM"]=0

  sol=cfba_moment(model,mod2,kcat,MW=mw,verbose=3,RHS=0.27,solver="glpkAPI",medval=3600*22.6) 

}% end dontrun
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FBA }
\keyword{ gene expression }
\keyword{ EFBA reaction }% __ONLY ONE__ keyword per line
