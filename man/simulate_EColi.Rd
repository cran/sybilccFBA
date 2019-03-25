\name{simulate_EColi}
\alias{simulate_EColi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Ecoli
}
\description{
calls different methods to simulate metabolism in EColi model.
}
\usage{
simulate_EColi(model, mod2, kcat, mw, budget_C, CSList, atpmz = TRUE, trns_rxns = NULL,
                cpx_stoich = NULL, solver = "glpkAPI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{   An object of class \code{\link{modelorg}}.}
  \item{mod2}{ An object of class \code{\link{modelorg}} with only irreversible reactions.
        It can be sent  to save time of recalculating it with each call.}
  \item{kcat}{ kcat values in unit 1/S. Contains three slots: reaction id,direction(dirxn),value(val)}
  \item{mw}{ list of molecular weights of all genes, using function calc_MW, in units g/mol}

  \item{budget_C}{  the budget C, for EColi 0.27}
  \item{CSList}{
    the list of carbon sources to be simulated, given as reaction id's
    }
  \item{atpmz}{
    set ATPM reaction to zero, default TRUE
    }
  \item{trns_rxns}{
    trnsport reactions to be excluded from budget (may contain other reactions if necessary)
    }
  \item{cpx_stoich}{
       giving the stoichiometry of complexes: data frame containing at least two columns 'genes','stoich' 
               'gene'  : delimited string of genes(ordered by geneid(e.g bnumber)),
               'stoich': number of subunits of each gene in the same order in 'genes'
    Default: NULL (e.g. all stoichiometry coefficients equal one)
   }

 \item{solver}{
       Single character string giving the solver package to use.  See
       \code{\link{SYBIL_SETTINGS}} for possible values.\cr
       Default: \code{SYBIL_SETTINGS("SOLVER")}.
     }
}
\details{
This function is to show examples of simulation
}
\value{
return a LIST, 
\item{gr_data}{the objective value (growth rate) of all carbon sources}
\item{uptake}{the uptake of the different conditions}
\item{all_flx}{data frame containing carbon source, reaction id, and flux in column flx}
\item{rg_MC}{reaction gene molecular crowding(MC), data frame containing carbon source, reaction id, 
  and flux in column flx, geneConc, rxnMC: MC for reaction, rgMC: MC for gene/complex,
  the result is sorted in decreasing order of MC}
\item{all_flx_MC}{reaction molecular crowding(MC), data frame containing carbon source, reaction id, 
  and flux in column flx, geneConc, rxnMC: MC for reaction, rgMC: MC for gene/complex,
  result not sorted.}
}
\references{
Desouki, Abdelmoneim. "Algorithms for improving the predictive power of flux balance analysis." PhD diss., 2016.
}
\author{
Abdelmoneim Amer Desouki
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{getccFBA_mat}},\code{\link{cfba_moment_mr}}
}
\examples{
 \dontrun{
    data(iML1515)
    #1-get metabolic model
    model=iML1515
    #2-get Molecular weights
    print(load(paste0(path.package("sybilccFBA"), '/extdata/mw_iML1515.RData')))

    mw=mw_iML1515
    mw=rbind(mw,data.frame(Synonym="s0001", mw=0.001))
    colnames(mw)[1]='gene'
    #3-get kcat list
    kl=read.csv(stringsAsFactors=FALSE,paste0(path.package("sybilccFBA"), 
        '/extdata/','allKcats_upd34_dd_h.csv'))
        kl=kl[!is.na(kl[,'ijo_id']),]
        kcat=data.frame(rxn_id=kl[,'ijo_id'],val=kl[,'kcat_max'],dirxn=kl[,'dirxn'],src=kl[,'src'],
        stringsAsFactors=FALSE)
        kcat=kcat[kcat[,'rxn_id']\%in\% react_id(model),]
        kcat[(is.na(kcat[,'src'])),'src']='Max'
    #4-get complex stoichiometry if used
    cpx_stoich =read.csv(paste0(path.package("sybilccFBA"), 
            '/extdata/','cpx_stoich_me.csv'),stringsAsFactors=FALSE)
    #5-identify Carbon sources to be tested
    csl=c("EX_glc__D(e)_b","EX_glyc(e)_b", "EX_ac(e)_b" ,  "EX_fru(e)_b",  
    "EX_pyr(e)_b",  "EX_gal(e)_b" ,
     "EX_lac__L(e)_b",  "EX_malt(e)_b" ,"EX_mal__L(e)_b", "EX_fum(e)_b" ,  
     "EX_xyl__D(e)_b","EX_man(e)_b", "EX_tre(e)_b",
     "EX_mnl(e)_b",   "EX_g6p(e)_b",   "EX_succ(e)_b", "EX_gam(e)_b",  
     "EX_sbt__D(e)_b", "EX_glcn(e)_b",
     "EX_rib__D(e)_b","EX_gsn(e)_b" ,"EX_ala__L(e)_b",  "EX_akg(e)_b" ,  
     "EX_acgam(e)_b")

    msrd=c(0.66,0.47,0.29,0.54,0.41,0.24,0.41,0.52,0.55,0.47,0.51,0.35,0.48,0.61,
        0.78,0.50,0.40,0.48,0.68,0.41,0.37,0.24,0.24,0.61)
    CA=c(6,3,2,6,3,6,3,12,4,4,5,6,12,6,6,4,6,6,6,5,10,3,5,8)

    #6-get irreversible model
    sum(react_id(model) \%in\% gsub('_b$','',csl))
    model1=model
    react_rev(model1)[react_id(model) \%in\% gsub('_b$','',csl)]=TRUE
    mod2=mod2irrev(model1)

     react_name(mod2)[react_id(mod2) \%in\% csl]
     
     uppbnd(mod2)[react_id(mod2) \%in\% csl]=0
     uppbnd(mod2)[react_id(mod2) \%in\% gsub('_b$','_f',csl)]=0

    uppbnd(mod2)[react_id(mod2)=="EX_o2(e)_b"]=1000
    trns_rxns=grepl("tex$",react_id(model))

    ##Call function
    tmp_res=simulate_EColi(model,mod2,mw=mw,budget_C=0.27,kcat=kcat,cpx_stoich=cpx_stoich,
    atpmz=FALSE,trns_rxns=trns_rxns,CSList=csl)
    bm=tmp_res[[1]]
    cor.test(bm[,'flx'],msrd,method='spearman')
     plot(msrd,bm[,'flx'],ylab="Predicted Growth Rate",xlim=c(0,0.8),
                  ylim=c(0,max(bm[,'flx'])),xlab='Measured Growth Rate',
        main=sprintf("Effect of Keff, Corr=\%.2f,nkcat=\%d",cor(bm[,'flx'],msrd),nrow(kcat)))
    abline(a=0,b=1,col="red",lwd=2,lty=2)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cost constraint FBA }% use one of  RShowDoc("KEYWORDS")
