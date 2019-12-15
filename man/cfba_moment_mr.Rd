\name{cfba_moment_mr}

\alias{cfba_moment_mr}

\encoding{utf8}

\title{ Function: cfba_moment_mr: implement MOMENT method 
                 taking into account multifunctional enzymes}
\description{
 This function uses GPR, kcat, and molecular weights to calculate fluxes 
 according to MOMENT method taking into account multifunctional enzymes. 
 Whenever a protein i was involved in more than one reaction, we introduced auxiliary
 concentration variables xi,j for each of these reactions. These xi,j replaced the 
 global concentration variable gi for the protein in the corresponding equation that 
 limits the flux through this reaction based on the enzyme concentration. 
 The sum of the xi,j is then equal to the total concentration of protein gi 
 included in the global enzyme solvent capacity constraint.
}

\usage{
cfba_moment_mr(model,mod2=NULL, Kcat,MW=NULL,
selected_rxns=NULL,verboseMode=2,objVal=NULL,
RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),C_mu_coef = 0,medval=NULL,
                 runFVA = FALSE, fvaRxn = NULL)
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
objective value(biomass)}

 \item{solver}{
       Single character string giving the solver package to use.  See
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("SOLVER")}.
     }
 \item{C_mu_coef}{ used to have C as a linear function of mu (biomass) : C = RHS + C_mu_coef*Biomass}
 \item{medval}{ median of Kcat values , used for missing values}
 \item{runFVA}{ flag to choose to run flux variability default FALSE}
 \item{fvaRxn}{ optional parameter to choose set of reaction ids to run FVA on them.
   Ids are from the irreversible model default all reactions. Ignored when runFVA 
   is not set.}

}
\details{
Main steps
    1- Add variables for all genes
    2- for each selected reaction: parse gpr,
    3- Add variables accordingly and constraints
	4- Add solvant constraint
}
\value{
    returns a list containing slots:
	 \item{sol}{solution of the problem.}
	 \item{prob}{object of class \code{\link{sysBiolAlg}} that contains the linear problem, 
			this can be used for further processing like adding more constraints. 
            To save it, function \code{\link{writeProb}} can be used.}
     \item{geneCol}{mapping of genes to variables in the problem.}
     \item{geneConc}{the concentration of each gene, when the gene is catalyzing more than
     one reaction there will be a row with 'rxn' column set to NA containing the total.}
     \item{rxnMC}{for each reaction (GPR) the molecular crowding of it (total sum to budget)}
      \item{rxnGeneMC}{the contribution of each gene to all of its reactions.}
 }

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{Adadi, R., Volkmer, B., Milo, R., Heinemann, M., & Shlomi, T. (2012). 
Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network 
with Kinetic Parameters, 8(7). doi:10.1371/journal.pcbi.1002575

Gelius-Dietrich, G., Desouki, A. A., Fritzemeier, C. J., & Lercher, M. J. (2013). 
sybil–Efficient constraint-based modelling in R. BMC systems biology, 7(1), 125.
}

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

  sol_mr=cfba_moment_mr(model,mod2,kcat,MW=mw,verbose=2,RHS=0.27,solver="glpkAPI",medval=3600*22.6) 
   bm_rxn = which(obj_coef(mod2)!=0)
   print(sprintf('biomass=\%f',sol_mr$sol$fluxes[bm_rxn]))
  # Enzyme concentrations:
     gconc=sol_mr$geneConc
  
}% end dontrun

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FBA }
\keyword{ MOMENT }
\keyword{ cost constraint FBA }% __ONLY ONE__ keyword per line
