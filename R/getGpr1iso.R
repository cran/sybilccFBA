## 22/10/2015

 getGpr1iso <- function (model,cpx_stoich,MW){
### TODO: add list of active genes/ know iso enzym/rxn pair
### preprocessing of model to get the minimal cost iso enzyme and eliminate the OR from GPR's
#   cpx_stoich: giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 1/10/2015
##               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
##               'stoich': number of subunits of each gene in the same order in 'genes'

#	MW measurement for gene using readfaa.r, calc_MW or gene annotation [units g/mmol]
        # Column 1 is the gene and column 2 is the MW
# Return: gpr1iso: new rules
#         rgst: matrix of reaction gene with stoichiometry
	
	nOR = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('|','',x,fixed=TRUE)))
	nAND = sapply(gprRules(model),function(x) nchar(x)-nchar(gsub('&','',x,fixed=TRUE)))

	### 0. Preprocessing evaluate complex MW's
	# mod_cmp=unlist(sapply(which(nOR==0 & nAND>0) ,function(x)paste(sort(allGenes(model)[rxnGeneMat(model)[x,]!=0]),collapse=',')))
	r1gflg=(rowSums(rxnGeneMat(model))==1)
	    mod_cpx=data.frame(stringsAsFactors=FALSE,rsn=which(r1gflg),rxn=react_id(model)[r1gflg],bool_expr=gpr(model)[r1gflg],
		 genes=gpr(model)[r1gflg],termNo=rep(1,sum(r1gflg)),nOR=nOR[r1gflg],nAND=nAND[r1gflg],gprtype=rep('r1g',sum(r1gflg)))
	
	r1cpxflg=(nOR==0 & nAND > 0)
	
	r1cpx_genes=sapply(which(r1cpxflg),function(x) return(paste(sort(allGenes(model)[rxnGeneMat(model)[x,]!=0]),collapse=',')))
	if(length(r1cpx_genes)>0){
		mod_cpx=rbind(mod_cpx,data.frame(stringsAsFactors=FALSE,rsn=which(r1cpxflg),rxn=react_id(model)[r1cpxflg],
					bool_expr=gpr(model)[r1cpxflg],
					genes=r1cpx_genes,termNo=1,nOR=nOR[r1cpxflg],nAND=nAND[r1cpxflg],gprtype='r1cpx'))
	}
	
	rg_ind = which(rxnGeneMat(model)!=0,arr.ind=T)
	rmgflg=(rg_ind[,1] %in% which(nOR > 0 & nAND == 0))
	if(any(rmgflg)){
		#termNo:used to choose when equal MW
		mod_cpx=rbind(mod_cpx,data.frame(stringsAsFactors=FALSE,rsn=rg_ind[rmgflg,1],rxn=react_id(model)[rg_ind[rmgflg,1]],
							bool_expr=gpr(model)[rg_ind[rmgflg,1]],
									genes=allGenes(model)[rg_ind[rmgflg,2]],termNo=1:sum(rmgflg),nOR=nOR[rg_ind[rmgflg,1]],nAND=nAND[rg_ind[rmgflg,1]],gprtype='rmg'))
	}
	## table(mod_cpx[,'gprtype'])
	rmcpxflg=(nOR > 0 & nAND > 0)

	mod_rmcpx=NULL
	xx= sapply(which(rmcpxflg),function(x) {
		rl=gpr(model)[x]
		rl=gsub("\\)"," ) ",rl)# 
		rl=gsub("\\("," ( ",rl)# 
		pr=unlist(strsplit(gsub(' OR ',' or ',rl,ignore.case=T)," or "))
		pr=gsub("[()]","",gsub(' AND ',' and ',pr,ignore.case=T))
		pr_genes=sapply(pr,function(z) paste(sort(gsub(' ','',unlist(strsplit(z,' and ')))),collapse=','))
		mod_rmcpx <<- rbind(mod_rmcpx,data.frame(stringsAsFactors=FALSE,rsn=x,rxn=react_id(model)[x],bool_expr=pr,genes=pr_genes,
					termNo=1:length(pr),nOR=nOR[x],nAND=nAND[x],gprtype='rmcpx'))
		return(length(pr))
	})

	mod_cpx =rbind(mod_cpx,mod_rmcpx)
    # browser()
	if(!is.null(cpx_stoich)){
		mod_cpx_st =data.frame(stringsAsFactors=FALSE,mod_cpx,
			stoich=cpx_stoich[match(mod_cpx[,'genes'],cpx_stoich[,'genes']),'stoich'])
	}else{
		mod_cpx_st =data.frame(stringsAsFactors=FALSE,mod_cpx,stoich=as.character(NA))
	}
	# browser();
	rownames(MW)=MW[,1]
	mod_cpx_MW1s=NULL # all ones
	mod_cpx_MWtxt=NULL
	mod_st_MW =NULL # using stoich.
	for(i in 1:nrow(mod_cpx_st)){
			genes=mod_cpx_st[i,'genes']
			stoich=mod_cpx_st[i,'stoich']
			mws=MW[unlist(strsplit(genes,',')),2]
			mod_cpx_MWtxt <- c(mod_cpx_MWtxt,paste(mws,collapse=','))
			mod_cpx_MW1s <- c(mod_cpx_MW1s,sum(mws))
			mod_st_MW <- c(mod_st_MW, sum(mws * as.numeric(unlist(strsplit(stoich,',')))))
	}

		termNG=sapply(mod_cpx_st[,'genes'],function(x) nchar(x)-nchar(gsub(',','',x,fixed=TRUE)))+1
	# gpr1iso= data.frame(stringsAsFactors=FALSE,gpr1iso,termNG)

	mod_cpx_mw = data.frame(stringsAsFactors=FALSE,mod_cpx_st,termNG,mod_cpx_MW1s,mod_st_MW,mod_cpx_MWtxt,
						MW=ifelse(is.na(mod_st_MW),mod_cpx_MW1s,mod_st_MW))
						
	if(!is.null(cpx_stoich))#rep(NA,nrow(mod_cpx_st)  # for debugging purposes
			mod_cpx_mw = data.frame(stringsAsFactors=FALSE,mod_cpx_mw,cpx_stoich[match(mod_cpx_mw[,'genes'],mod_cpx_mw[,'genes']),]);#c('enzMW','MW1s','cpx_MWtxt')

### Preprocessing 2. choose best iso enzyme for each rxn  --------------------------------				
	flg=mod_cpx_mw[,'nOR']>0
	rxn_iso_cpx=by(mod_cpx_mw[mod_cpx_mw[,'nOR']>0,'MW'],mod_cpx_mw[mod_cpx_mw[,'nOR']>0,'rsn'],min)
	## each iso gpr : choose term with min(MW)
	gpr1iso= mod_cpx_mw[mod_cpx_mw[,'nOR']==0 | paste(mod_cpx_mw[,'rsn'],mod_cpx_mw[,'MW']) %in% paste(names(rxn_iso_cpx),as.vector(rxn_iso_cpx)),]
	#13/3/2017: problem of equal iso enzymes
	if(any(duplicated(gpr1iso[,'rsn']))){
		equal_iso=by(gpr1iso[gpr1iso[,'nOR']>0,'termNo'],gpr1iso[gpr1iso[,'nOR']>0,'rsn'],min)
		gpr1iso= gpr1iso[gpr1iso[,'nOR']==0 | paste(gpr1iso[,'rsn'],gpr1iso[,'termNo']) %in% 
		                         paste(names(equal_iso),as.vector(equal_iso)),]
	}
	# termNG=sapply(gpr1iso[,'genes'],function(x) nchar(x)-nchar(gsub(',','',x,fixed=TRUE)))+1
	# gpr1iso= data.frame(stringsAsFactors=FALSE,gpr1iso,termNG)
	gpr1iso[is.na(gpr1iso[,'stoich']) & gpr1iso[,'termNG']==1,'stoich']=1
	rgst=matrix(0,ncol=ncol(rxnGeneMat(model)),nrow=nrow(rxnGeneMat(model)))
	rownames(rgst)=react_id(model)
	colnames(rgst)=allGenes(model)
	rgst[cbind(gpr1iso[gpr1iso[,'termNG']==1,'rxn'],gpr1iso[gpr1iso[,'termNG']==1,'genes'])] = gpr1iso[gpr1iso[,'termNG']==1,'stoich']
	for(i in which(gpr1iso[,'termNG']>1)){
		if(is.na(gpr1iso[i,'stoich'])){
			rgst[gpr1iso[i,'rxn'],unlist(strsplit(gpr1iso[i,'genes'],','))]= 1
            # browser()
            gpr1iso[i,'stoich']=paste(rep(1,gpr1iso[i,'termNG']),collapse=',')
		}else{
			rgst[gpr1iso[i,'rxn'],unlist(strsplit(gpr1iso[i,'genes'],','))]=as.numeric(unlist(strsplit(gpr1iso[i,'stoich'],',')))
		}
	}
	
	###
	return(list(gpr1iso=gpr1iso,mod_cpx_mw=mod_cpx_mw,rgst=rgst))
}

