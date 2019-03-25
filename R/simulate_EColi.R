# 2/12/2018

simulate_EColi<-function(model,mod2,kcat,mw,budget_C,CSList,atpmz=TRUE,trns_rxns=NULL,
             cpx_stoich=NULL,solver='glpkAPI'){
#trns_rxns: trnsport reactions to be excluded
#CSList: the list of carbon sources
# atpmz: set ATPM reaction to zero

if(atpmz){
    lowbnd(mod2)[react_id(mod2)=="ATPM"]=0
    uppbnd(mod2)[react_id(mod2)=="ATPM"]=0
}    
 nr=react_num(mod2)
 medianVal=median(kcat[,'val'])
 mod2R=mod2
if(!is.null(trns_rxns)){ 
    selected_rxns=react_id(model)[gpr(model)!="" & !trns_rxns]
    print(sprintf('excluded transport reactions:%d',sum(trns_rxns)))
}else{
    selected_rxns=react_id(model)[gpr(model)!=""]
}

    mrMat = getccFBA_mat(model,mod2,kcat,MW=mw,verboseMode=2,RHS=budget_C,cpx_stoich=cpx_stoich,medval=3600*medianVal,selected_rxns=selected_rxns)
    rgm=mrMat$rxnGeneCol
    cnames=mrMat$cnames
    gpr1iso=mrMat$gpr1iso
 # bm_rxn=which(obj_coef(mod2)!=0)
    solverParm=NA
    # gr_data = NULL
    all_flx=NULL
    all_flx_MC=NULL
    all_rg_MC=NULL
  for(cs in 1:length(CSList)){
		print(paste0(cs,' ',CSList[cs]))
		mod2=mod2R
		uppbnd(mod2)[react_id(mod2) %in% CSList]=0
		uppbnd(mod2)[react_id(mod2)==CSList[cs]]=1000
	
	prob=sysBiolAlg(mod2,LHS=mrMat$LHS,rlb=mrMat$rlb,rub=unlist(mrMat$rub),rtype=mrMat$rtype,lb=mrMat$clb,ub=mrMat$cub,obj=mrMat$obj_cf,
			lpdir='max',cnames=mrMat$cnames,solver=solver,  algorithm="ccFBA",solverParm=solverParm,
					pname=sprintf("ccFBA mfe %s: %s,cpxst,kcatupd25",solver,mod_name(model)),
					writeProbToFileName=sprintf('EC_ccFBA_Mat_%s.lp',solver));
	sol=optimizeProb(prob);
	str(sol)
	print(getMeanStatus(sol$stat,solver))
    ###
    	all_flx=rbind(all_flx,data.frame(stringsAsFactors=F,cs,csname=CSList[cs],rxn_id=react_id(mod2),flx=sol$fluxes[1:nr],kcat=mrMat$kcat_ir,kcat_src=mrMat$kcat_src_ir,
					enzConc=sol$fluxes[mrMat$rxnExprCol[match(1:nr,mrMat$rxnExprCol[,'flxCol']),'exprCol']],obj=obj_coef(mod2),ub=uppbnd(mod2),ubR=uppbnd(mod2R)))
    flx_MC = data.frame(stringsAsFactors=FALSE,rgm, #[,c('rsn','rxn_id','gene','gMW','gpr.MW')],
      flx=sol$fluxes[rgm[,'rxnCol']],
      grConc=ifelse(is.na(rgm[,'grCol']),0,sol$fluxes[rgm[,'grCol']]),
	enzConc=sol$fluxes[rgm[,'gpr.exprCol']]/as.numeric(rgm[,'stoich']),#shouldn't be NA # stoich shouldn't be 0:in case of r1g
	geneConc=ifelse(is.na(rgm[,'gCol']),0,sol$fluxes[rgm[,'gCol']])
	)

	flx_MC=cbind(flx_MC,rgMC=rgm[,'gMW'] * flx_MC[,'grConc'] ,#/ as.numeric(flx_MC[,'stoich'])
			rxnMC=rgm[,'gpr.MW'] * flx_MC[,'enzConc'] ) 
	# rg_MC=flx_MC[flx_MC[,'rgMC']>1e-7 & flx_MC[,'flx']>1e-7 ,]
	rg_MC=flx_MC#25/5/2017
	rg_MC=rg_MC[order(rg_MC[,'rgMC'],decreasing=T),]
	rg_MC=cbind(rgMCorder=1:nrow(rg_MC),rg_MC)
	rg_MC=rg_MC[order(rg_MC[,'rxnMC'],decreasing=T),]
	rg_MC=cbind(rxnMCorder=1:nrow(rg_MC),rg_MC)
	# rgMCorder=order(as.numeric(rg_MC[,'rgMC']),decreasing=T)
	# rxnMCorder=order(as.numeric(rg_MC[,'rxnMC']),decreasing=T)
	all_flx_MC=rbind(all_flx_MC,cbind(cs,csname=CSList[cs],flx_MC))
	all_rg_MC=rbind(all_rg_MC,cbind(cs,csname=CSList[cs],rg_MC))

}

    upt=all_flx[all_flx[,'csname']==all_flx[,'rxn_id'],]
    gr_data=all_flx[all_flx[,'obj']!=0,]

    return(list(gr_data=gr_data,uptake=upt,all_flx=all_flx,rg_MC=all_rg_MC,all_flx_MC=all_flx_MC))
}