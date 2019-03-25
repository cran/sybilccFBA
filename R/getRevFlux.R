getRevFlux <- function(model,irrxns,fdirrev){
# for reversible rxns get (_f-_b) / the members of irrevmod class can be used instead!1/6/2015
fluxes=NULL;
 for (r in(react_id(model))){
	if(! react_rev(model)[react_id(model)==r]){
	 fluxes=rbind(fluxes,data.frame(rxn=r,fwd=fdirrev[which(irrxns==r)],bwd=0))
	}else{
		fluxes=rbind(fluxes,data.frame(rxn=r,fwd=fdirrev[which(irrxns==paste(r,'_f',sep=""))]
	               ,bwd=fdirrev[which(irrxns==paste(r,'_b',sep=""))]))
	
	} 
}
return(fluxes)
}
