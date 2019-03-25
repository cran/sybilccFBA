#given  Kcat vector (in 1/Sec), and molecular weights  (Maximize biomass)
#1-Add gi vars (gene concentrations)
#2-add gpr constraints: add f/b for reversible rxns react_rev
#		vi/kcati <= gi
#       vi cat by g1 & g2 :  add aux variable then Aux <= all(gi), vi<= Kcat*Aux
#       vi cat by g1 | g2 :   OR is eliminated by gpr1iso, choose the smallest enzyme
#      complex gpr recursively
##... ----------
#3-	mr:multiple reactions (multiple function enzymes
# 		for genes catalyzing more n rxn(n>1): add gij: j 1:n  such that Sum(gij) <= gi
#4-add Density constraint Sum(MW*gi)<=C

getccFBA_mat <- function (model, mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,objVal=NULL,
									RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL,
							
		cpx_stoich = NULL,C_mu_coef=0){
#Multiple OR are not done pairwise but evaluated in one line in contrast to MOMENT original implementation.
#addCols
#	Kcat measurement for set of rxns  (Unit: 1/s)
#	MW measurement for gene using readfaa.r, calc_MW or gene annotation
#   RHS : C  g protein/gDW  (metabolic)
#   objVal: when set the minimum C that gives this growth, will be returned as objective
#   cpx_stoich: giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 1/10/2015
##               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
##               'stoich': number of subunits of each gene in the same order in 'genes'
##  C_mu_coef: used to have C as a linear function of mu (biomass) : C = RHS + C_mu_coef*Biomass

#			units: mmol/gDW/h
	# vi : mmol/gDW/hr  , rxn catal by gi with Kcati
	# Kcat: hr-1
	# gi: mmol/gDW
	# MWi: g/mmol (val of AA in g/mol, divide by 1000)
	# Calc_MW
	# C : g/gDW protein (metabolic)
	# Linear effect on flux, can be optimized to min error 0.27

	Kcat[,'val']<-as.numeric(Kcat[,'val'])*60*60 ;# convert to 1/hr
	MW[,2]=as.numeric(MW[,2]); #g/mol instead of g/mmol
	kcat_ir = NULL  #9/8/2015
	if(length(mod2)==0){ ## 28/10/2015
			mod2=mod2irrev(model)
		}
	if(length(medval)==0) medval=median(Kcat[,'val'])

	if(length(selected_rxns)==0)
		{selected_rxns=react_id(model)[gpr(model)!="" ]} #& gpr(model)!="s0001" exclude s0001 in 25 rxns +4 rxns

	# nOR = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('|','',x,fixed=TRUE)))
	# nAND = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('&','',x,fixed=TRUE)))

	### 0. Preprocessing evaluate complex MW's
	# mod_cmp=unlist(sapply(which(nOR==0 & nAND>0) ,function(x)paste(sort(allGenes(model)[rxnGeneMat(model)[x,]!=0]),collapse=',')))
	
	print("Preprocessing isoenzymes ... ");
	r1i=getGpr1iso(model,MW=MW,cpx_stoich=cpx_stoich)
	gpr1iso=r1i$gpr1iso;
	mod_cpx_mw=r1i$mod_cpx_mw; rgst=r1i$rgst;
    if(any(is.na(gpr1iso[,'stoich']))){
        stop("There are NA values in gpr1iso")
    }
### start making constraint matrix ---------------------------------------------------
	print("Preparing LP ... ");
	LHS <- S(mod2)
	n=react_num(mod2)
	m=met_num(mod2)
	rnames=paste('S',met_id(mod2),sep='_')
	rlb=rub=rep(0,m);
	rtype=rep('E',m);
	cnames=paste('R',react_id(mod2),sep='_')
	colid=n+1 #index of next col
	rowind=m+1
	aux_id=1
	geneCol=NULL;#column in problem corresponding to given gene
	rxnExprCol = NULL
	
	
# add all variables once
#Add variables for all genes : geneCol
	rgst[! react_id(model)%in% selected_rxns,]=0
	selectedGenes = (colSums(rgst!=0)!=0)
	ng = sum( selectedGenes )
	geneind=rep(NA,ncol(rgst))  ## use NA for unused genes and not 0[ returns different lengths]
    geneind[selectedGenes]= ncol(LHS) + (1 :  ng)
	names(geneind) = allGenes(model)
	genecnt = colSums(rgst!=0);
	rg_ind = which(rgst!=0,arr.ind=T)
	geneCol=data.frame(stringsAsFactors=FALSE,gene=allGenes(model)[selectedGenes],rxn=NA,Col = geneind[selectedGenes],cnt=genecnt[selectedGenes],
					vname=paste(ifelse(genecnt[selectedGenes]==1,'g1','gmr'),allGenes(model)[selectedGenes], genecnt[selectedGenes],sep='_'), r=NA,g=which(selectedGenes))# geneind[selectedGenes],
								
	geneCol[geneCol[,'cnt']==1,'r']= rg_ind[match(geneCol[geneCol[,'cnt']==1,'g'],rg_ind[,2]),1]
    # LHS <- cBind(LHS,matrix(0,ncol=nrow(geneCol),nrow=nrow(LHS)))
    LHS <- cbind(LHS,matrix(0,ncol=nrow(geneCol),nrow=nrow(LHS)))
	cnames=c(cnames,geneCol[,'vname'])
	# browser()
	rg_mfe = rg_ind[rg_ind[,2] %in% which(colSums(rgst!=0) > 1) ,]
	if(nrow(rg_mfe)>0){
		mfeCol = data.frame(stringsAsFactors=FALSE, gene=allGenes(model)[rg_mfe[,2]], rxn=react_id(model)[rg_mfe[,1]],
			 Col= ncol(LHS) + (1:nrow(rg_mfe)), cnt=genecnt[rg_mfe[,2]], rsn=rg_mfe[,1],g=rg_mfe[,2])
			
		mfeCol= data.frame(stringsAsFactors=FALSE,mfeCol, vname = paste('gr',mfeCol[,'gene'],mfeCol[,'rxn'],sep='_'))#,mfeCol[,'Col']
		# LHS <- cBind(LHS,matrix(0,ncol=nrow(mfeCol),nrow=nrow(LHS)))
		LHS <- cbind(LHS,matrix(0,ncol=nrow(mfeCol),nrow=nrow(LHS)))
		cnames = c(cnames,mfeCol[,'vname'])
	}else{mfeCol=NULL}
	
	mgrflg=which(rowSums(rgst!=0) > 1 & react_id(model) %in% selected_rxns) # Add auxiliary variables to multiple gene rules
	if(length(mgrflg)>0){
	        gprCol = data.frame(stringsAsFactors=FALSE,       # containing more than one gene
			       rsn=mgrflg,rxn=react_id(model)[mgrflg],Col= ncol(LHS) + (1:length(mgrflg)),
				   vname = paste('cpx',react_id(model)[mgrflg],sep='_')#,ncol(LHS) + (1:length(mgrflg))
			)
	}else{
		gprCol=NULL
	}
	
	if(!is.null(gprCol)){
		# LHS <- cBind(LHS,matrix(0,ncol=nrow(gprCol),nrow=nrow(LHS)))
		LHS <- cbind(LHS,matrix(0,ncol=nrow(gprCol),nrow=nrow(LHS)))
		cnames = c(cnames,gprCol[,'vname'])
	}

    exprCol=ifelse(gpr1iso[,'rsn'] %in% gprCol[,'rsn'],gprCol[match(gpr1iso[,'rsn'], gprCol[,'rsn']),'Col'],# r1cpx
         ifelse(gpr1iso[,'genes'] %in% mfeCol[,'gene'], mfeCol[match(gpr1iso[,'rsn'], mfeCol[,'rsn']),'Col'], #r1g:mfe
                 geneCol[match(gpr1iso[,'genes'],geneCol[,'gene']),'Col']))									 #r1g:1to1
	   table(substring(cnames[exprCol],1,2)) 
		gpr1iso=data.frame(stringsAsFactors=FALSE,gpr1iso,exprCol=exprCol)
	# }else{
		# gpr1iso=data.frame(stringsAsFactors=FALSE,gpr1iso,exprCol=NA)
		# browser()
	# }
	# if(!is.null(exprCol))
	### ----------------------------------------
	# vind,vind_b,kval,rkval
	kcat_pv=kcat_nv=rep(medval,react_num(model))
	# browser()
	kcat_pv[match(Kcat[Kcat[,"dirxn"]==1 & Kcat[,'rxn_id'] %in% react_id(model),"rxn_id"],react_id(model))]=
	              Kcat[Kcat[,"dirxn"]==1 & Kcat[,'rxn_id'] %in% react_id(model),"val"]
	kcat_nv[match(Kcat[Kcat[,"dirxn"]==-1 & Kcat[,'rxn_id'] %in% react_id(model),"rxn_id"],react_id(model))]=
	                    Kcat[Kcat[,"dirxn"]==-1 & Kcat[,'rxn_id'] %in% react_id(model),"val"]

	kcat_src_pv=kcat_src_nv=rep("medianVal",react_num(model))
	kcat_src_pv[match(Kcat[Kcat[,"dirxn"]==1 & Kcat[,'rxn_id'] %in% react_id(model),"rxn_id"],react_id(model))]=
	        Kcat[Kcat[,"dirxn"]==1 & Kcat[,'rxn_id'] %in% react_id(model),"src"]
	kcat_src_nv[match(Kcat[Kcat[,"dirxn"]==-1 & Kcat[,'rxn_id'] %in% react_id(model),"rxn_id"],react_id(model))]=
	Kcat[Kcat[,"dirxn"]==-1 & Kcat[,'rxn_id'] %in% react_id(model),"src"]

	kcat_ir=rep(medval,react_num(mod2))
	kcat_ir[grepl("_b$",react_id(mod2))]=kcat_nv[irrev2rev(mod2)[grepl("_b$",react_id(mod2))]]
	kcat_ir[!grepl("_b$",react_id(mod2))]=kcat_pv[irrev2rev(mod2)[!grepl("_b$",react_id(mod2))]]# not backward
	# kcat_ir=kcat_ir[match(selected_rxns,react_id(mod2))] #10/3/2017
	
	kcat_src_ir=rep("medianVal",react_num(mod2))
	kcat_src_ir[grepl("_b$",react_id(mod2))]=kcat_src_nv[irrev2rev(mod2)[grepl("_b$",react_id(mod2))]]
	kcat_src_ir[!grepl("_b$",react_id(mod2))]=kcat_src_pv[irrev2rev(mod2)[!grepl("_b$",react_id(mod2))]]
	
	### Add constraint
	#AddRow Vind - Kcat*colind <= 0 {all selected rxns}
	### ----------------------------------------
	# selected_rxns =selected_rxns[!selected_rxns %in% gpr1iso[gpr1iso[,'genes']=='s0001','rxn']]#grepl('s0001',gpr1iso[,'genes'])
	# browser()
	rxns_ind = gpr1iso[match(selected_rxns,gpr1iso[,'rxn']),'rsn'] ##
	rxnkcatRow = data.frame(stringsAsFactors=FALSE,rsn=rxns_ind, rxn_id = selected_rxns, 
					 vind=ifelse(!react_rev(model)[rxns_ind],match(selected_rxns,react_id(mod2)),
					 rev2irrev(mod2)[rxns_ind,1]),
	                 vind_b=ifelse(!react_rev(model)[rxns_ind],NA,rev2irrev(mod2)[rxns_ind,2]),
					 gpr1iso[match(rxns_ind,gpr1iso[,'rsn']),c('stoich','termNG','exprCol')])
	rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow,kval=kcat_ir[rxnkcatRow[,'vind']],
					r1gstoich=as.numeric(ifelse(rxnkcatRow[,'termNG']>1,1,rxnkcatRow[,'stoich'])))
	if(any(!is.na(rxnkcatRow[,'vind_b']))){ #11/3/2017
		rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow, rkval=kcat_ir[rxnkcatRow[,'vind_b']])
		}else{
			rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow, rkval=NA)
		}
	rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow,rowind=nrow(LHS)+(1:length(rxns_ind)) )
	## 1-rkcat_ constr
	# browser()
	consNames=paste('rkcat',react_id(mod2)[rxnkcatRow[,'vind']],sep='_')
	# LHS <- rBind(LHS,matrix(0,nrow=nrow(rxnkcatRow),ncol=ncol(LHS)))
	LHS <- rbind(LHS,matrix(0,nrow=nrow(rxnkcatRow),ncol=ncol(LHS)))
	LHS[as.matrix(rxnkcatRow[,c('rowind','vind')])]=1
	LHS[as.matrix(rxnkcatRow[,c('rowind','exprCol')])]= -1*rxnkcatRow[,'kval']/rxnkcatRow[,'r1gstoich']
	rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow,consNames)
	rnames=c(rnames,consNames)
	rlb=c(rlb,rep(0,nrow(rxnkcatRow))); rub=c(rub,rep(0,nrow(rxnkcatRow)));
	rtype=c(rtype,rep('U',nrow(rxnkcatRow)))
	###
	
	flg=!is.na(rxnkcatRow[,'vind_b'])
	rowind_b=consNames_b=rep(NA,nrow(rxnkcatRow))
	if(sum(flg)>0){
		consNames_b[flg]=paste('rbkcat',react_id(mod2)[rxnkcatRow[flg,'vind_b']],sep='_')
		rowind_b[flg]=nrow(LHS)+(1:sum(flg))
		rxnkcatRow = data.frame(stringsAsFactors=FALSE, rxnkcatRow,consNames_b,rowind_b)
		# LHS <- rBind(LHS,matrix(0,nrow=sum(flg),ncol=ncol(LHS)))
		LHS <- rbind(LHS,matrix(0,nrow=sum(flg),ncol=ncol(LHS)))
		LHS[as.matrix(rxnkcatRow[flg,c('rowind_b','vind_b')])]=1
		LHS[as.matrix(rxnkcatRow[flg,c('rowind_b','exprCol')])]= -1*rxnkcatRow[flg,'rkval']/rxnkcatRow[flg,'r1gstoich']
		rnames=c(rnames,consNames_b[flg])
		rlb=c(rlb,rep(0,sum(flg))); rub=c(rub,rep(0,sum(flg)));
		rtype=c(rtype,rep('U',sum(flg)))
	}
	
	rxnExprCol=rbind(cbind(flxCol=rxnkcatRow[,'vind'],exprCol=rxnkcatRow[,'exprCol'],rowind=rxnkcatRow[,'rowind']),
	                 cbind(flxCol=rxnkcatRow[flg,'vind_b'],exprCol=rxnkcatRow[flg,'exprCol'],rowind=rxnkcatRow[flg,'rowind_b']))
	###------------------------------------------------------------------------
	## 2-min_ constr for complexes termNG > 1, for each gene: mapping and to minimum
	if(verboseMode>=2) print('Adding constraints for mapping AND to min');
	# browser();
	 rxns_ind=which(rowSums(rgst!=0)>1 & (react_id(model) %in% selected_rxns))
	rg_r1cpx = rg_ind[rg_ind[,1] %in%  rxns_ind,]
	minGiRow = data.frame(stringsAsFactors=FALSE,rsn=rg_r1cpx[,1],gsn=rg_r1cpx[,2],rxn=react_id(model)[rg_r1cpx[,1]],
					gene=allGenes(model)[rg_r1cpx[,2]],stoich=as.numeric(rgst[rg_r1cpx]),exprCol=gpr1iso[match(rg_r1cpx[,1],gpr1iso[,'rsn']),'exprCol'])
	minGiRow = data.frame(stringsAsFactors=FALSE,minGiRow,
	       geneCol=ifelse(minGiRow[,'gene'] %in% mfeCol[,'gene'], 
					mfeCol[match(paste(minGiRow[,'rsn'],minGiRow[,'gene']), paste(mfeCol[,'rsn'],mfeCol[,'gene'])),'Col'], 
								geneCol[match(minGiRow[,'gene'],geneCol[,'gene']),'Col']) )
	if(nrow(minGiRow)>0){
		rowind=nrow(LHS)+(1:nrow(minGiRow))
		consNames=paste('min',minGiRow[,'gene'],minGiRow[,'rxn'],minGiRow[,'stoich'],sep='_')
		minGiRow = data.frame(stringsAsFactors=FALSE,minGiRow,rowind,consNames)
		# LHS <- rBind(LHS,matrix(0,nrow=nrow(minGiRow),ncol=ncol(LHS)))
		LHS <- rbind(LHS,matrix(0,nrow=nrow(minGiRow),ncol=ncol(LHS)))
		LHS[cbind(minGiRow[,'rowind'],minGiRow[,'exprCol'])] = -1*minGiRow[,'stoich']
		LHS[cbind(minGiRow[,'rowind'],minGiRow[,'geneCol'])] = 1
		rnames=c(rnames,consNames)
		rlb=c(rlb,rep(0,nrow(minGiRow))); rub=c(rub,rep(0,nrow(minGiRow)));
		# rtype=c(rtype,rep('E',nrow(minGiRow)))
		rtype=c(rtype,rep('L',nrow(minGiRow)))#should be U, 29/3/2017
	}	
	## 3-mfe_ constr---------------------------------------------------------------------
	#Add multiple reactions constraints
if(verboseMode>=2) print('Adding multiple reactions constraints');
# browser()
	mfeRow = allGenes(model)[colSums(rgst!=0)>1]
	if(length(mfeRow)>0){
		rowind=nrow(LHS) + (1:length(mfeRow))
		mfeRow =data.frame(stringsAsFactors=FALSE,gene=mfeRow,rowind,Col=geneCol[match(mfeRow,geneCol[,'gene']),'Col'])
		mfeCol=data.frame(stringsAsFactors=FALSE,mfeCol,grow=mfeRow[match(mfeCol[,'gene'],mfeRow[,'gene']),])
		# LHS <- rBind(LHS,matrix(0,nrow=nrow(mfeRow),ncol=ncol(LHS)))
		LHS <- rbind(LHS,matrix(0,nrow=nrow(mfeRow),ncol=ncol(LHS)))
		LHS[as.matrix(mfeCol[,c('grow.rowind','Col')])] = 1  ## geneRxn
		LHS[as.matrix(mfeRow[,c('rowind','Col')])] = -1  ## total gene
		consNames=paste('mfe',mfeRow[,'gene'],sep='_')#,mfeRow[,'rowind']
		rnames=c(rnames,consNames)
			rlb=c(rlb,rep(0,nrow(mfeRow))); rub=c(rub,rep(0,nrow(mfeRow)));
			rtype=c(rtype,rep('E',nrow(mfeRow)))
	}		
##############II. Add MW constraint------------------
## Sum(gi. MWi)<= C[gdwpr/gDW]
	fluxind = c(1:n)
	# geneind = as.numeric(geneCol[match(allGenes(model),geneCol[is.na(geneCol[,"rxn"]),'gene']),"Col"]) #
	clb=cub=rep(0,ncol(LHS))
	clb[fluxind] =lowbnd(mod2)
	cub[fluxind] =uppbnd(mod2)
	cub[(n+1):ncol(LHS)] = 1000;  # C is less than 1, can be 2017

if(length(MW )>0 ){
	if(verboseMode>=2) print('Adding capacity constraint.');
	# browser()
	gMW = MW[match(allGenes(model), MW[,1]),2]
	if(sum((!is.na(geneind)) & is.na(gMW)) > 0) 
			print(sprintf('Warning : some genes have no MW:%d',sum((!is.na(geneind)) & is.na(gMW)) ));
	lcind = geneind[(!is.na(geneind)) & (!is.na(gMW)) & geneind!=0]
	lnz = gMW[(!is.na(geneind)) & (!is.na(gMW)) & geneind!=0]
	
	if(!is.null(RHS) ){	#test value of upper bound (C)
	# Add solvant constraint	
		# LHS=rBind(LHS,0)
		LHS=rbind(LHS,0)
		rowind=nrow(LHS)
		LHS[rowind,lcind] = lnz
		LHS[rowind,fluxind[obj_coef(mod2)!=0]] = -C_mu_coef ## 3/10/2015 - 22/9/2016		
		rlb=c(rlb,0);      rub=c(rub,RHS);
		rnames=c(rnames,'cc')
		rtype=c(rtype,'U')
		cc_ind = rowind

		obj_cf =rep(0,ncol(LHS))
		obj_cf[fluxind]=obj_coef(mod2)
		mod_objc_ind =NA;		
		obj_dir = 'max';# should be as model objective
		
	} else {  # Add constraint on FBA objective
		if(is.null(objVal) ){	#calculate minimum cost and gene concentration
			tmp_sol = optimizeProb(mod2,solver=solver,retOptSol=FALSE);
			if(tmp_sol$ok!=0 ||  length(checkSolStat(stat=tmp_sol$stat,solver=solver))!=0 ) {#check if solution is feasible
					if(is.na(checkSolStat(stat=tmp_sol$stat,solver=solver)[1])) print("couldn't check solution status");
					stop(sprintf("Failed to calculate objective value, sol status: %s", getMeanStatus(code=tmp_sol$stat,solver=solver) ));
			}	
		
			objVal = tmp_sol$obj;
			print(sprintf('Objective value to be used: %f',objVal))
		}
		
		if(is.na(objVal)){ 
				stop('objVal parameter is NA.')
		}
		# objc=getObjCoefs(problem(prob),j=c(1:n))
		objc=obj_coef(mod2)#getObjCoefs(problem(prob),j=c(1:n))
		# cind=which(objc!=0)
		cind=fluxind[objc!=0]
		nzval=objc[cind]
		# LHS=rBind(LHS,0)
		LHS=rbind(LHS,0)

		LHS[rowind,cind] = nzval
		rlb=c(rlb,objVal);      rub=c(rub,Inf);
		rnames=c(rnames,'ObjC_')
		rtype=c(rtype,'U')
		mod_objc_ind=rowind;
		cc_ind = NA;
		obj_dir = 'min';# minimize
    	 rowind=rowind+1;  
		# changeObjCoefs(lp = problem(prob),j=c(1:n),obj_coef=rep(0,n))
		obj_cf=rep(0,ncol(LHS))
		obj_cf[lcind] = lnz;
	}
	}
			 #######################################################
	rg_ind2 = which(rxnGeneMat(mod2)!=0,arr.ind=TRUE)  # 121,12m1,12
	rxnGeneCol=data.frame(stringsAsFactors=FALSE,rsn=rg_ind2[,1],gsn=rg_ind2[,2],
										rxn_id=react_id(mod2)[rg_ind2[,1]],  gene=allGenes(mod2)[rg_ind2[,2]],
										revRxn=irrev2rev(mod2)[rg_ind2[,1]],
										rNG=rowSums(rgst!=0)[irrev2rev(mod2)[rg_ind2[,1]]], genercnt=genecnt[rg_ind2[,2]],
										rxnCol=fluxind[rg_ind2[,1]], gCol=geneind[rg_ind2[,2]],  
										gpr=gpr1iso[match(irrev2rev(mod2)[rg_ind2[,1]],gpr1iso[,'rsn']),c('exprCol','MW','mod_cpx_MW1s','genes')],
										kcat=kcat_ir[rg_ind2[,1]],kcat_src=kcat_src_ir[rg_ind2[,1]],
										gMW=gMW[rg_ind2[,2]], stoich=rgst[cbind(irrev2rev(mod2)[rg_ind2[,1]],rg_ind2[,2])]
										)
											
			rxnGeneCol=data.frame(stringsAsFactors=FALSE,rxnGeneCol,
						grCol=ifelse(genecnt[rg_ind2[,2]] > 1, 
							mfeCol[match(paste(irrev2rev(mod2)[rg_ind2[,1]],rg_ind2[,2]), paste(mfeCol[,'rsn'],mfeCol[,'g'])),'Col'], 
							geneind[rg_ind2[,2]]),
							exprCname=cnames[rxnGeneCol[,'gpr.exprCol']],gColname=cnames[rxnGeneCol[,'gCol']]
							)
							
			rxnGeneCol=data.frame(stringsAsFactors=FALSE,rxnGeneCol,grColname=cnames[rxnGeneCol[,'grCol']])
						
###
## cnames: rxn_id,g_gene,cmplx_rxn,g_gene_rxn,iso_rxn
	out <- list( 	LHS = LHS,
					rlb = rlb, rub=rub, rnames=rnames, rtype=rtype,   	# rows
					clb = clb, cub=cub, cnames=cnames,     				# columns
					fluxind = fluxind, geneind = geneind,
					rxnkcatRow=rxnkcatRow,minGiRow=minGiRow,mfeRow=mfeRow,
					rgst=rgst,mod_cpx_mw=mod_cpx_mw,                      # complex stoichiometry
					rxnExprCol = rxnExprCol,       			              # flxCol,exprCol,rowind,reaction total expression
					rxnGeneCol = rxnGeneCol,       			              # rxn_id,gene,Col : used to get MCrowding
					geneCol = geneCol, mfeCol=mfeCol,gpr1iso=gpr1iso,     # low level structure, used for debugging
					kcat_ir = kcat_ir, kcat_src_ir=kcat_src_ir, geneMW = gMW,
					model = mod2,
					obj_cf = obj_cf,    obj_dir =obj_dir,
					mod_objc_ind = mod_objc_ind,            # row index of constraint on FBA model objective
					cc_ind = cc_ind                         # capacity constraint index
					# C   # budget rub[cc_ind]
				);
				
	return(out);
}
