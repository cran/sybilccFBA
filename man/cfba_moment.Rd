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
RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL,
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
	 \item{sol}{solution of the problem, instance of class \code{\link{optObj}}.}
	 \item{prob}{object of class \code{\link{sysBiolAlg}} that contains the linear problem, 
			this can be used for further processing like adding more constraints. 
            To save it, function \code{\link{writeProb}} can be used.}
     \item{geneCol}{mapping of genes to variables in the problem.}
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

  sol=cfba_moment(model,mod2,kcat,MW=mw,verbose=2,RHS=0.27,solver="glpkAPI",medval=3600*22.6) 
  bm_rxn = which(obj_coef(mod2)!=0)
  print(sprintf('biomass=\%f',sol$sol$fluxes[bm_rxn]))
 # Enzyme concentrations:
 
}% end dontrun
data(Ec_core)
model=Ec_core
genedef=read.csv(paste0(path.package("sybilccFBA"), '/extdata/Ec_core_genedef.csv'),
     stringsAsFactors=FALSE)
mw=data.frame(gene=genedef[,'gene'],mw=genedef[,'mw'],stringsAsFactors=FALSE)
mw[mw[,1]=='s0001','mw']=0.001#spontenious
##########
##Kcats
kl=read.csv(stringsAsFactors=FALSE,paste0(path.package("sybilccFBA"), 
'/extdata/','allKcats_upd34_dd_h.csv'))
    kl=kl[!is.na(kl[,'ijo_id']),]
    kcat=data.frame(rxn_id=kl[,'ijo_id'],val=kl[,'kcat_max'],dirxn=kl[,'dirxn'],
                 src=kl[,'src'],stringsAsFactors=FALSE)
    kcat=kcat[kcat[,'rxn_id']\%in\% react_id(model),]
    kcat[(is.na(kcat[,'src'])),'src']='Max'
 ########## ----------------
 
mod2=mod2irrev(model)
uppbnd(mod2)[react_id(mod2)=="EX_o2(e)_b"]=1000
lowbnd(mod2)[react_id(mod2)=="ATPM"]=0
uppbnd(mod2)[react_id(mod2)=="ATPM"]=0 

nr=react_num(mod2)
medianVal=median(kcat[,'val'])
# sum(is.na(genedef[,'mw']))

C_mr=0.13#as not all genes exist
sim_name=paste0('Ec_org_med',round(medianVal,2),'_C',100*C_mr)

cpx_stoich =read.csv(paste0(path.package("sybilccFBA"), 
'/extdata/','cpx_stoich_me.csv'),stringsAsFactors=FALSE)


##-----------------
CSList=c("R_EX_glc_e__b","R_EX_glyc_e__b","R_EX_ac_e__b","R_EX_fru_e__b",
"R_EX_pyr_e__b","R_EX_gal_e__b",
"R_EX_lac_L_e__b","R_EX_malt_e__b","R_EX_mal_L_e__b","R_EX_fum_e__b","R_EX_xyl_D_e__b",
"R_EX_man_e__b","R_EX_tre_e__b",
"R_EX_mnl_e__b","R_EX_g6p_e__b","R_EX_succ_e__b","R_EX_gam_e__b","R_EX_sbt_D_e__b",
"R_EX_glcn_e__b",
"R_EX_rib_D_e__b","R_EX_gsn_e__b","R_EX_ala_L_e__b","R_EX_akg_e__b","R_EX_acgam_e__b")

msrd=c(0.66,0.47,0.29,0.54,0.41,0.24,0.41,0.52,0.55,0.47,0.51,0.35,0.48,0.61,
0.78,0.50,0.40,0.48,0.68,0.41,0.37,0.24,0.24,0.61)
CA=c(6,3,2,6,3,6,3,12,4,4,5,6,12,6,6,4,6,6,6,5,10,3,5,8)

CSList=substring(gsub('_e_','(e)',CSList),3)
 react_name(mod2)[react_id(mod2) \%in\% CSList]
 msrd=msrd[CSList \%in\% react_id(mod2) ]    
 CA=CA[CSList \%in\% react_id(mod2)]  
 CSList=CSList[CSList \%in\% react_id(mod2)]
 uppbnd(mod2)[react_id(mod2) \%in\% CSList]=0
 mod2R=mod2
##---------------------
bm_rxn=which(obj_coef(mod2)!=0)
all_flx=NULL
all_flx_MC=NULL
all_rg_MC=NULL
Kcatorg=kcat
solver='glpkAPI'
solverParm=NA

for(cs in 1:length(CSList)){
		print(CSList[cs])
		mod2=mod2R
		uppbnd(mod2)[react_id(mod2) \%in\% CSList]=0
		uppbnd(mod2)[react_id(mod2)==CSList[cs]]=1000
	sol_org=cfba_moment(model,mod2,kcat,MW=mw,verbose=2,RHS=0.27,solver="glpkAPI",
    medval=3600*medianVal) 
  
	### preparing output -------------------
	all_flx=rbind(all_flx,data.frame(stringsAsFactors=FALSE,cs,csname=CSList[cs],
            rxn_id=react_id(mod2),
            flx=sol_org$sol$fluxes[1:nr], ub=uppbnd(mod2),ubR=uppbnd(mod2R)))
        # print(paste0("nrow all_rg_MC=",nrow(all_rg_MC)))
}

upt=all_flx[all_flx[,'csname']==all_flx[,'rxn_id'],]
bm=all_flx[react_id(mod2)[obj_coef(mod2)!=0]==all_flx[,'rxn_id'],]

cor.test(bm[,'flx'],msrd,method='spearman')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FBA }
\keyword{ MOMENT }
\keyword{ cost constraint FBA }% __ONLY ONE__ keyword per line
