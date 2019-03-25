## Abdelmoneim Amer Desouki, 12/5/2015



#------------------------------------------------------------------------------#
#                 definition of the class sysBiolAlg_ccFBA                      #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_ccFBA",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_ccFBA
setMethod(f = "initialize",
          signature = "sysBiolAlg_ccFBA",
          definition = function(.Object,
								model,
                                LHS,
                                rlb,rub,rtype,
								lb,ub,obj,#ctype,
								lpdir,
								cnames=NULL,
								rnames=NULL,
								pname = NULL,
                                scaling = NULL,
                                writeProbToFileName = NULL, ...) {

						stopifnot(nrow(LHS)==length(rlb),
									nrow(LHS)==length(rub),
									nrow(LHS)==length(rtype),
									ncol(LHS)==length(lb),
									ncol(LHS)==length(ub),
									# ncol(LHS)== length(ctype),
									ncol(LHS)== length(obj)
									);
									
									## model should be modelorg and Irrev is removed
									lb[1:react_num(model)]=lowbnd(model)
									ub[1:react_num(model)]=uppbnd(model)
									# use the default constructor for class sysBiolAlg
                  
				  .Object <- callNextMethod(.Object,
                                            sbalg      = "ccFBA",
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = c(1:ncol(LHS)),
                                            nCols      = ncol(LHS),
                                            nRows      = nrow(LHS),
                                            mat        = LHS,
                                            ub         = ub,
                                            lb         = lb,
                                            obj        = obj,
                                            rlb        = rlb,
                                            rtype      = rtype,
                                            lpdir      = lpdir,
                                            rub        = rub,
                                            # ctype      = ctype,
                                            cnames     = cnames,
                                            rnames     = rnames,
                                            pname     = pname,
                                            algPar     = list(),
                                            ...)

                  # write problem to lp file
                  if (!is.null(writeProbToFileName)) {
                      writeProb(problem(.Object),
                                fname = as.character(writeProbToFileName))
                  }

                  validObject(.Object)
                  return(.Object)
          }
)
