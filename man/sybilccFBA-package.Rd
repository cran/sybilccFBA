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
\examples{
\dontrun{
data(Ec_core)
model=Ec_core
genedef=read.csv(paste0(path.package("sybilccFBA"), '/extdata/Ec_core_genedef.csv'),
     stringsAsFactors=FALSE)
mw=data.frame(gene=genedef[,'gene'],mw=genedef[,'mw'],stringsAsFactors=FALSE)
mw[mw[,1]=='s0001','mw']=0.001#spontenious

kl=read.csv(stringsAsFactors=FALSE,paste0(path.package("sybilccFBA"), 
            '/extdata/','allKcats_upd34_dd_h.csv'))
    kl=kl[!is.na(kl[,'ijo_id']),]
    kcat=data.frame(rxn_id=kl[,'ijo_id'],val=kl[,'kcat_max'],dirxn=kl[,'dirxn'],src=kl[,'src'],
                stringsAsFactors=FALSE)
    kcat=kcat[kcat[,'rxn_id']\%in\% react_id(model),]
    kcat[(is.na(kcat[,'src'])),'src']='Max'

mod2=mod2irrev(model)
uppbnd(mod2)[react_id(mod2)=="EX_o2(e)_b"]=1000
lowbnd(mod2)[react_id(mod2)=="ATPM"]=0
uppbnd(mod2)[react_id(mod2)=="ATPM"]=0 


nr=react_num(mod2)
medianVal=median(kcat[,'val'])
# sum(is.na(genedef[,'mw']))

C=0.13#as not all genes exist
sim_name=paste0('Ec_core_med',round(medianVal,2),'_C',100*C)

cpx_stoich =read.csv(paste0(path.package("sybilccFBA"), 
                '/extdata/','cpx_stoich_me.csv'),stringsAsFactors=FALSE)

CSList=c("R_EX_glc_e__b","R_EX_glyc_e__b","R_EX_ac_e__b","R_EX_fru_e__b",
                "R_EX_pyr_e__b","R_EX_gal_e__b",
    "R_EX_lac_L_e__b","R_EX_malt_e__b","R_EX_mal_L_e__b","R_EX_fum_e__b",
                "R_EX_xyl_D_e__b","R_EX_man_e__b","R_EX_tre_e__b",
    "R_EX_mnl_e__b","R_EX_g6p_e__b","R_EX_succ_e__b","R_EX_gam_e__b",
                "R_EX_sbt_D_e__b","R_EX_glcn_e__b",
    "R_EX_rib_D_e__b","R_EX_gsn_e__b","R_EX_ala_L_e__b",
    "R_EX_akg_e__b","R_EX_acgam_e__b")

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
# xt=FALSE
st=TRUE
selected_rxns=react_id(model)[gpr(model)!=""]

if(st){
    mrMat = getccFBA_mat(model,mod2,kcat,MW=mw,verbose=2,RHS=C,cpx_stoich=cpx_stoich,
                              medval=3600*medianVal,selected_rxns=selected_rxns)
}else{#nost
    mrMat = getccFBA_mat(model,mod2,kcat,MW=mw,verbose=2,RHS=C,cpx_stoich=NULL,
                          medval=3600*medianVal,selected_rxns=selected_rxns)
 }
mrMat_st = mrMat

r1i=getGpr1iso(model,MW=mw,cpx_stoich=cpx_stoich)
	gpr1iso=r1i$gpr1iso;
	mod_cpx_mw=r1i$mod_cpx_mw; rgst=r1i$rgst;
  sum(is.na(gpr1iso[,'stoich']))

nOR = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('|','',x,fixed=TRUE)))
nAND = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('&','',x,fixed=TRUE)))

  rgm=mrMat$rxnGeneCol
write.csv(file='rgm_aa.csv',cbind(rgm,nOR[rgm[,'revRxn']],gpr(model)[rgm[,'revRxn']]))

rgm=rgm[!is.na(rgm[,'gCol']) & rgm[,'stoich']!=0,]# only used genes
cnames=mrMat$cnames
gpr1iso=mrMat$gpr1iso
bm_rxn=which(obj_coef(mod2)!=0)
all_flx=NULL
all_flx_mr=NULL
all_flx_org=NULL

Kcatorg=kcat
solver='glpkAPI'
solverParm=NA

for(cs in 1:length(CSList)){
		print(CSList[cs])
		mod2=mod2R
		uppbnd(mod2)[react_id(mod2) \%in\% CSList]=0
		uppbnd(mod2)[react_id(mod2)==CSList[cs]]=1000
	
	prob=sysBiolAlg(mod2,LHS=mrMat$LHS,rlb=mrMat$rlb,rub=mrMat$rub,rtype=mrMat$rtype,
       lb=mrMat$clb,ub=mrMat$cub,obj=mrMat$obj_cf,lpdir='max',cnames=mrMat$cnames,solver=solver,
               algorithm="ccFBA",solverParm=solverParm,
			pname=sprintf("ccFBA mfe \%s: \%s,cpxst,kcatupd25",solver,mod_name(model)),
					writeProbToFileName=sprintf('EC_ccFBA_Mat_\%s.lp',solver));
	sol=optimizeProb(prob);
	str(sol)
	print(getMeanStatus(sol$stat,solver))
    sol_mr=cfba_moment_mr(model,mod2,kcat,MW=mw,verbose=1,RHS=C,solver="glpkAPI",
             medval=3600*medianVal)
	sol_org=cfba_moment(model,mod2,kcat,MW=mw,verbose=1,RHS=C,solver="glpkAPI",
             medval=3600*medianVal)
    ### preparing output -------------------
	all_flx=rbind(all_flx,data.frame(stringsAsFactors=FALSE,cs,csname=CSList[cs],
                                rxn_id=react_id(mod2),flx=sol$fluxes[1:nr]))
    all_flx_mr=rbind(all_flx_mr,data.frame(stringsAsFactors=FALSE,cs,csname=CSList[cs],
                        rxn_id=react_id(mod2),flx=sol_mr$sol$fluxes[1:nr]))
    all_flx_org=rbind(all_flx_org,data.frame(stringsAsFactors=FALSE,cs,
            csname=CSList[cs],rxn_id=react_id(mod2),flx=sol_org$sol$fluxes[1:nr]))
}

upt=all_flx[all_flx[,'csname']==all_flx[,'rxn_id'],]

#bm=all_flx[all_flx[,'obj']!=0,]
bm=all_flx[react_id(mod2)[obj_coef(mod2)!=0]==all_flx[,'rxn_id'],]
bm_mr=all_flx_mr[react_id(mod2)[obj_coef(mod2)!=0]==all_flx_mr[,'rxn_id'],]
bm_org=all_flx_org[react_id(mod2)[obj_coef(mod2)!=0]==all_flx_org[,'rxn_id'],]

cor.test(bm[,'flx'],msrd,method='spearman')
cor.test(bm_mr[,'flx'],msrd,method='spearman')
cor.test(bm_org[,'flx'],msrd,method='spearman')

bm=cbind(bm,msrd)

#Wong
cor(as.numeric(bm[,"flx"]),(as.numeric(bm[,"flx"])/as.numeric(upt[,"flx"]))^0.5 )

\dontrun{
 
 plot(msrd,bm[,'flx'],ylab="pred bm",xlim=c(0,0.8),ylim=c(0,max(bm_org[,'flx'])),
   pch=19,col='red',
                    main="Predicting growth rate using Ec_core")
abline(a=0,b=1,col="grey",lwd=2,lty=2)
points(msrd,bm_org[,'flx'],ylab="pred bm",xlim=c(0,0.8),ylim=c(0,max(bm_org[,'flx'])),
   pch=19,col='green')
points(msrd,bm_mr[,'flx'],ylab="pred bm",xlim=c(0,0.8),ylim=c(0,max(bm_org[,'flx'])),
   pch=19,col='darkblue')

legend(legend=c('MOMENT','MOMENT_mr','ccFBA st'),col=c('green','darkblue','red'),
                x=0,y=1,cex=0.75,pch=19)

################################top 5 genes####################
# selected_genes =c('b2323','b2472','b3774','b2779','b0431') # excl 'b0241', too many rxns 248
selected_genes =c('b1773','b0118','b1779') # excl 'b0241', too many rxns 248
	sel_rxns =react_id(mod2)[(rowSums(rxnGeneMat(mod2)[,allGenes(mod2)\%in\% selected_genes])>0)]
	sel_rxns=c(sel_rxns,CSList[cs],react_id(mod2)[obj_coef(mod2)!=0])

	
 write.csv(file=sprintf('iAF_top4g_flx_le\%s.csv',sim_name),
                    all_flx[all_flx[,'rxn_id']\%in\% sel_rxns,])
}
}
}
\keyword{ package }
\seealso{
	\code{\link{sybil}} 
	\code{\link{cfba_moment}} 
}