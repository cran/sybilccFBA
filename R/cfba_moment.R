#given  Kcat vector (in 1/Sec), and molecular weights
#1-Add gi vars
#2-add gpr constraints: add f/b for reversible rxns react_rev??
#3-add Density constraint
#... ----------
cfba_moment <- function (model,mod2=NULL, Kcat,MW=NULL,selected_rxns=NULL,verboseMode=2,objVal=NULL,
		RHS=NULL,solver=SYBIL_SETTINGS("SOLVER"),medval=NULL,runFVA=FALSE,fvaRxn=NULL){
#Multiple OR are not done pairwise but evaluated in one line.
#addCols
#	Kcat measurement for set of rxns  (unit: 1/s)
#	MW measurement for gene using readfaa.r
#require(sybil)
geneCol=NULL;#column in problemcoresponding to given gene

Kcat[,2]<-as.numeric(Kcat[,2])*60*60 ;# convert to 1/hr
MW[,2]=as.numeric(MW[,2]); #g/mol instead of g/mmol
if(length(mod2)==0){
		mod2=mod2irrev(model)
	}
if(length(medval)==0) medval=median(Kcat[,2])

if(length(selected_rxns)==0)
    {selected_rxns=react_id(model)[gpr(model)!="" & gpr(model)!="s0001"]} #exclude s0001 in 25 rxns +4 rxns

 print("Preparing LP ... ");
prob <- sysBiolAlg(mod2, algorithm = "fba",solver=solver)
n=react_num(mod2)
m=met_num(mod2)

colid=n+1
rowind=m+1
aux_id=1
# add all variables once
#Add variables for all genes
	for( g in allGenes(mod2)){
		if(g!="s0001"){
					geneCol=rbind(geneCol,cbind(gene=g,Col=colid))
					addCols(lp = problem(prob),1)
					changeColsBnds(lp = problem(prob),colid,lb=0,ub=1000)
					colid=colid+1;			
				}}

#ng=length(allGenes(model))
rxnkcatRow=NULL;
for (r in selected_rxns){
	rl=gpr(model)[react_id(model)==r]
	#rl=gsub("s0001",)
	# remove s0001 gene from rules in iAF model
	if (rl=='( b0875  or  s0001 )') rl='b0875';
	if (rl=='( b1377  or  b0929  or  b2215  or  b0241  or  s0001  or  b3875  or  b1319  or  b0957 )')
		 rl='( b1377  or  b0929  or  b2215  or  b0241 or  b3875  or  b1319  or  b0957 )';
	if (rl=='( s0001  or  b0451 )') rl='(b0451)';
	if (rl=='( s0001  or  b3927 )') rl='( b3927 )';
	#ijo model
	if(rl=="(b3927 or s0001)") rl="b3927";
	if(rl=="(s0001 or b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)" )
	     rl="(b0957 or b3875 or b2215 or b0241 or b1319 or b1377 or b0929)"
    if(rl=="(s0001 or b0875)") rl="b0875";
	if(rl=="(b0451 or s0001)") rl="b0451";
		#dh  additions for iJO
	if(rl=="(b1377 or b0241 or b3875 or s0001 or b0957 or b2215 or b1319 or b0929)"){rl<-"(b1377 or b0241 or b3875 or b0957 or b2215 or b1319 or b0929)"}
	#dh  additions for iML
	if(rl=="(b1377 or b2215 or b1319 or b3875 or s0001 or b0957 or b0241 or b0929)") rl="(b1377 or b2215 or b1319 or b3875 or b0957 or b0241 or b0929)";
	if(rl=="(s0001 or b0451)") rl="b0451";
	if(rl=="(s0001 or b3927)") rl="b3927";
	if(rl=="(s0001 or b1779)") rl="b1779";
	if(rl=="(b1779 or s0001)") rl="b1779";
    #Ec_core
	if(rl=="(b0875 or s0001)") rl="b0875";
    
	if (verboseMode > 2){ print(r)}
	if(r %in% Kcat[Kcat[,"dirxn"]==1,1]){
		kval=as.numeric(Kcat[Kcat[,"dirxn"]==1 & Kcat[,1]==r,"val"]);
		}else{kval=medval;}
	
	if(r %in% Kcat[Kcat[,"dirxn"]==-1,1]){
		rkval=as.numeric(Kcat[Kcat[,"dirxn"]==-1 & Kcat[,1]==r,"val"]);
		}else{rkval=medval;}
	
	if(!react_rev(model)[react_id(model)==r]){ #if ! reversible
			vind=which(react_id(mod2)==r)
			vind_b=0
			lnz=1
		}else{# Add row with median value to rev direction
			#vind=c(which(react_id(mod2)==paste(r,"_f",sep="")),which(react_id(mod2)==paste(r,"_b",sep="")))		
			#lnz=c(1,1)
			vind=which(react_id(mod2)==paste(r,"_f",sep=""))
			vind_b=which(react_id(mod2)==paste(r,"_b",sep=""))		
			lnz=1
		}
	rl=gsub("\\)"," ) ",rl)# 
	rl=gsub("\\("," ( ",rl)# 
	rl=gsub(" OR "," or ",rl)#21/3/2017
	rl=gsub(" AND "," and ",rl)
		# to be sure that 'and' will be a whole word and not part of any gene name
	pr=lapply(strsplit(unlist(strsplit(rl," or "))," and "),function(x) gsub("[() ]","",x))
	if (verboseMode > 2){
		print(pr)
		print(length(pr[[1]]))
	}
	if( length(pr)==1) {# no OR (only one term) 
	#add row vi<=gi*kcat
		if(length(pr[[1]])==1){
#########################################1-one Gene----------------
		    gene=pr[[1]]
		#get the corresponding col of the gene
		colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
		if (verboseMode > 2){print(list(vind,colind))}
	
		#AddRow Vind - Kcat*colind <=0
		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind,colind)),
              nzval = list(c(lnz,-1*kval))
			  ,rnames = r
			 )
			 rowind=rowind+1;
			 ## kcatRow: store rows of kcat values  25/1/2017
             # if(r=='ACALDtpp')browser()
			 rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind, rxn_id = r, dirxn=1,
							rowind=rowind-1,colind,kval));
	        if(vind_b>0){
		  		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind_b,colind)),
              nzval = list(c(lnz,-1*rkval))
			  ,rnames = paste(r,"_b",sep="")
			 )
			 rowind=rowind+1;            
			 rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind_b, rxn_id = r, dirxn=-1,
							rowind=rowind-1,colind,kval=rkval));
		  }			 
		}else{#straight ANDs
#########################################2-straight ANDs----------------
		if(verboseMode>2) print(sprintf("Rxn:%s gpr: %s, straight ANDs",r,rl));
		
		#add aux variable
		addCols(lp = problem(prob),1)
		aux_id=colid
		colid=colid+1;			
		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind,aux_id)),
              nzval = list(c(lnz,-1*kval))
			  ,rnames = r
			 )
			 rowind=rowind+1;
			#31/1/2017			 
			 rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind, rxn_id = r, dirxn=1,
							rowind=rowind-1,colind=aux_id,kval));
	        		
		if(react_rev(model)[react_id(model)==r]){
    		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(vind_b,aux_id)),
              nzval = list(c(lnz,-1*rkval))
			  ,rnames = paste(r,"_b",sep="")
			 )
			 rowind=rowind+1;            
			 #31/1/2017
			rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind_b, rxn_id = r, dirxn=-1,
							rowind=rowind-1,colind=aux_id,kval=rkval));
		}
		for(gene in pr[[1]]){			
				colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
			addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(aux_id,colind)),
              nzval = list(c(1,-1))
			#  ,rnames = paster
			 )
			 rowind=rowind+1;   
		}
		}
		}else{#OR/AND
#########################################3-straight OR----------------
		#add row for the rxn
		   row_vals=rep(0,colid+length(pr))
			row_vals[vind]=lnz
			row_vals_b=rep(0,colid+length(pr))
			row_vals_b[vind_b]=lnz
		 for( p in 1:length(pr)){#ORed gene list
           if( length(pr[[p]])==1 ){
				gene=pr[[p]]
				colind=as.numeric(geneCol[geneCol[,"gene"]==gene,"Col"])
          	    #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=colind,val=-1*kval)
				row_vals[colind]=-1*kval
			#31/1/2017
     		    rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind, rxn_id = r, dirxn=1,
							rowind=NA,colind=colind,kval));
				if(react_rev(model)[react_id(model)==r]){
				 #changeMatrixRow(lp=problem(prob),i=rxn_row_id_b,j=colind,val=-1*rkval)
				  row_vals_b[colind]=-1*rkval
				  #31/1/2017
				  rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind_b, rxn_id = r, dirxn=-1,
									rowind=NA,colind=colind,kval=rkval));
				}
		  }else{#AND
		  #add auxiliary var then add rows
#########################################4-mixed OR and AND (sum of products)----------------		  
		if(verboseMode > 2) print(sprintf("Rxn:%s gpr: %s, auxiliary ANDs",r,rl));
			addCols(lp = problem(prob),1)
			  aux_id=colid
			colid=colid+1;			
          
			#changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=aux_id,val=-1*kval)	
			row_vals[aux_id]=-1*kval
			 #31/1/2017
			 # browser()
     		 rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind, rxn_id = r, dirxn=1,
							rowind=NA,colind=aux_id,kval));
			if(react_rev(model)[react_id(model)==r]){			
			   # changeMatrixRow(lp=problem(prob),i=rxn_row_id_b,j=aux_id,val=-1*rkval)
				row_vals_b[aux_id]=-1*rkval			   
								#31/1/2017
				rxnkcatRow = rbind(rxnkcatRow,data.frame(stringsAsFactors=FALSE,rsn=vind_b, rxn_id = r, dirxn=-1,
									rowind=NA,colind=aux_id,kval=rkval));
			}
			for(g in pr[[p]]){
                colind=as.numeric(geneCol[geneCol[,"gene"]==g,"Col"])
        	
		     addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(c(aux_id,colind)),          nzval = list(c(1,-1))
			  ,rnames = paste("Aux_",aux_id,sep="")
			 )
						 rowind=rowind+1;  
			}
          if (verboseMode > 2){ print(p)}
		  }		  
		  #changeMatrixRow(lp=problem(prob),i=rxn_row_id,j=which(row_vals!=0),val=row_vals[which(row_vals!=0)])	
		}  
			addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
            cind = list(which(row_vals!=0)),        nzval = list(row_vals[which(row_vals!=0)])
			    ,rnames = r
			 )
			#rxn_row_id=rowind
			rowind=rowind+1;
			rxnkcatRow[is.na(rxnkcatRow[,'rowind']) & rxnkcatRow[,'rxn_id']==r & rxnkcatRow[,'dirxn']==1,'rowind']=rowind-1;
			
			if(react_rev(model)[react_id(model)==r]){
             		addRowsToProb(lp = problem(prob),
              i = rowind ,              type = "U",
              lb = 0,              ub = 0,
              cind = list(which(row_vals_b!=0)),        nzval = list(row_vals_b[which(row_vals_b!=0)])
			  ,rnames = paste(r,"_b",sep="")
			 )
			#rxn_row_id_b=rowind
			rowind=rowind+1;  
			rxnkcatRow[is.na(rxnkcatRow[,'rowind']) & rxnkcatRow[,'rxn_id']==r & rxnkcatRow[,'dirxn']==-1,'rowind']=rowind-1;
	       }
		
				 #rowind=rowind+1;            
	}
}	
############change column bounds ##########
#n +1 ,colid
#changeColsBnds(lp = problem(prob),c(n+1:colid-1),lb=rep(0,colid-n),ub=rep(1000,colid-n))
###############II. Add MW constraint------------------
## Sum(gi * MWi) <= C[gdwpr/gDW]

if(length(MW )>0 ){
   lnz=NULL
   lcind=NULL
   # print(colid)
   if (verboseMode > 2){	print('Setting column bounds ...') }
		
   for(v in c((n+1):(colid-1))){
   # print(v);
		changeColsBnds(lp = problem(prob),v,lb=0,ub=1000)
	}
	
   for(g in geneCol[,"gene"]) {
		colind=as.numeric(geneCol[geneCol[,"gene"]==g,"Col"])
		#changeColsBnds(lp = problem(prob),colind,lb=0,ub=1000)
		if(g %in% MW[,1]){# consider only genes with computed MW.
			lnz=cbind(lnz,MW[MW[,1]==g,2])
			if (verboseMode > 2){print(c(g,MW[MW[,1]==g,2]))}
			lcind=cbind(lcind,colind)
			}
		}


	if(!is.null(RHS) ){	#test value of upper bound (C)
	# Add solvant constraint	
	   addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "U",
				  lb = 0,              ub = RHS,#C denotes the total weight of proteins,
				  cind = list(lcind),          nzval = list(lnz)
				  ,rnames = "MW"
				 )
	} else {
		if(is.null(objVal) ){	#calculate minimum cost and gene concentration
			tmp_sol = optimizeProb(mod2,solver=solver,retOptSol=FALSE);
			if(tmp_sol$ok!=0 ||  length(checkSolStat(stat=tmp_sol$stat,solver=solver))!=0 ) {#check if solution is feasible
					if(is.na(checkSolStat(stat=tmp_sol$stat,solver=solver)[1])) print("couldn't check solution status");
					stop(sprintf("Failed to calculate objective value, sol status: %s", getMeanStatus(code=tmp_sol$stat,solver=solver) ));
			}	
		
			objVal = tmp_sol$obj;
			print(sprintf('Objective value to be used: %f',objVal))
		}
		objc=getObjCoefs(problem(prob),j=c(1:n))
		cind=which(objc!=0)
		nzval=objc[cind]
		addRowsToProb(lp = problem(prob),
				  i = rowind ,              type = "L",
				  lb = objVal,              ub = Inf, # ub ignored in type L
				  cind = list(cind),          nzval = list(nzval)
				  ,rnames = paste("ObjC_",rowind,sep="")
				 )
							 rowind=rowind+1;  
		changeObjCoefs(lp = problem(prob),j=c(1:n),obj_coef=rep(0,n))
		#print("old obj..")
		#set new obj function
		changeObjCoefs(lp = problem(prob),j=as.numeric(lcind),obj_coef=lnz)
		setObjDir(lp = problem(prob),"min")
	}
	}
			 #######################################################
	if (verboseMode > 2) {                      
				            fname=format(Sys.time(), "moment_%Y%m%d_%H%M.lp");
					       print(sprintf("Writing the problem to: %s/%s...",getwd(),fname));
		                	writeProb(lp=problem(prob), fname)
							if(solver=="cplexAPI"){
								Nnz=NULL
								nr=getNumRows(problem(prob))
								lp=problem(prob)
								for(r in c(1:nr)){
								  # print(r)
									 rv=cplexAPI::getRowsCPLEX(env = lp@oobj@env, lp = lp@oobj@lp, begin = r-1, end = r-1)
									 #write rv,matind,val
									 for(cl in 1:length(rv$matind)){
										   Nnz=rbind(Nnz,cbind(i=r,j=rv$matind[cl],val=rv$matval[cl]))
										}
								}
								write.csv(file="nnz.csv",Nnz)
							}
                }
	############################
print("Solving ... ");	
	sol=optimizeProb(prob)
	
	## Check solution status  18/6/2015
	if(sol$ok!=0 ||  length(checkSolStat(stat=sol$stat,solver=solver))!=0 ) {#check if solution is feasible
		    if(is.na(checkSolStat(stat=sol$stat,solver=solver)[1])) print("couldn't check solution status");
			warning(sprintf("Solution is not optimal: %s", getMeanStatus(code=sol$stat,solver=solver) ));
	}
	
	actual_C = matrix(lnz,nrow=1) %*% matrix(sol$fluxes[lcind],ncol=1);
	if(!runFVA){
		return(list(sol=sol,geneCol=geneCol,prob=prob,actual_C=actual_C,rxnkcatRow=rxnkcatRow))
	## running FVA (6/8/2015)
	### ------------------- FVA ----------------------
	}else{
	if (verboseMode > 1) {print("Calculating FVA ...")}
	main_flx = sol$fluxes[1:react_num(mod2)]
	tol = SYBIL_SETTINGS("TOLERANCE");
	# upt_rxn=which(react_id(mod2)=='R_EX_ac_e__b')		
	bm_rxn = which(obj_coef(mod2)==1)
	## set biomass
	main_sol = sol
	# main_flx = sol$fluxes[1:react_num(mod2)]
	model_objVal = main_flx[bm_rxn] #sol$obj; #Check if not calculating C
	#------
	print(sprintf('Main sol obj rxn: %s, objVal:%f',react_id(mod2)[bm_rxn],main_flx[bm_rxn]))
	
	# objc=getObjCoefs(problem(prob),j=c(1:n))
	# cind=which(objc!=0)

	changeColsBnds(lp = problem(prob), bm_rxn, lb=model_objVal-tol, ub=model_objVal+tol) ## allow no deviation (it can only decrease)
	changeObjCoefs(lp = problem(prob),j=bm_rxn,obj_coef=0)
	
	bLP = backupProb(lp = problem(prob))
	
	 # if (verboseMode > 1) { progr <- sybil::.progressBar() }
	# intReact = c(100:115)#c(sample(which(main_flx>0.1 & gpr(mod2)!="" ),15),sample(which(main_flx < 1e-6 & gpr(mod2)!=""),15))
	  if(!is.null(fvaRxn)){intReact = fvaRxn}
	  else{ intReact = 1:react_num(mod2)}
	  
	  
		nObj <- 2 * length(intReact) #ifelse(isTRUE(modIrr), length(intReact), 2 * length(intReact)) 

        obj   <- numeric(nObj)
        ok    <- integer(nObj)
        stat  <- integer(nObj)
		
		oppositFlx    <- numeric(nObj)
		bmFlx    	  <- numeric(nObj)  # stores biomass value at min/max (effective when not fixed)
    
    for (i in seq(along = intReact)) {
		# the opposite rxn should carry no flux: the LP cost will insure this
			
		prob1=prob
		prob@problem@oobj = bLP@oobj
		bLP <- backupProb(bLP)
		remove(prob1)
			
			solpl <-  i	
			sol <- optimizeProb(prob, lpdir = "min",
                                react = intReact[i], obj_coef = 1)

				obj[solpl]   <- sol$obj
                ok[solpl]    <- sol$ok
                stat[solpl]  <- sol$stat
				bmFlx[solpl]  <- sol$fluxes[bm_rxn]
				
				if (matchrev(mod2)[intReact[i]] !=0 ){
						oppositFlx[solpl]    <- sol$fluxes[matchrev(mod2)[intReact[i]]];
				}
			
			if(sol$stat != 5){## TOBE removed
					print('bad status')
					# browser()
					}
			
			# ----------------------------------
			
			solpl <-  length(intReact) + i
        
            sol <- optimizeProb(prob, lpdir = "max",
                                react = intReact[i], obj_coef = 1)
                obj[solpl]   <- sol$obj
                ok[solpl]    <- sol$ok
                stat[solpl]  <- sol$stat
				bmFlx[solpl]  <- sol$fluxes[bm_rxn]
				
				if (matchrev(mod2)[intReact[i]] !=0 ){
						oppositFlx[solpl]    <- sol$fluxes[matchrev(mod2)[intReact[i]]];
				}
					
				if(sol$stat != 5){
					print('bad status')
					# browser()
					}
			#---------------------------
			
				# print(sprintf('rxn:%s,flx:%f, min:%f,max:%f,oposit:(%f,%f), biomass:(%f,%f)',react_id(mod2)[intReact[i]],main_flx[intReact[i]],
				# obj[solpl/2] ,obj[solpl],
				# oppositFlx[solpl/2],oppositFlx[solpl],bmFlx[solpl/2],bmFlx[solpl]))
        # if (verboseMode > 1) {
            # progr <- sybil::.progressBar(i, length(intReact), progr)
        # }

    }
	
	fva <- cbind (ok =ok, stat=stat, lpdir=c(rep('min',length(intReact)),rep('max',length(intReact))),fv=obj, mainFlx=main_flx[c(intReact,intReact)],
							oppositFlx=oppositFlx, bmFlx=bmFlx, rxns=c(intReact,intReact), rxnid=react_id(mod2)[c(intReact,intReact)],
						gpr=gpr(mod2)[c(intReact,intReact)],oposid=matchrev(mod2)[c(intReact,intReact)]);
	
	return(list(sol=main_sol,geneCol=geneCol,prob=prob,actual_C=actual_C,fva=fva));
}
}