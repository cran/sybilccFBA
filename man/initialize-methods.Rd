\name{initialize-methods}

\docType{methods}
\encoding{utf8}

\alias{initialize,sysBiolAlg_ccFBA-method}


%\alias{sysBiolAlg_ccFBA,sysBiolAlg_ccFBA-method}

\title{Initialize Problem Object}

\description{
  Initialize ccFBA Problem.
}

\usage{
\S4method{initialize}{sysBiolAlg_ccFBA}(.Object, model, LHS, rlb, rub, rtype, lb, ub, obj,
                 lpdir, cnames = NULL, rnames = NULL, pname = NULL,
                 scaling = NULL, writeProbToFileName = NULL, ...)

}

%\section{Methods}{
 % \describe{
  \arguments{
 \item{.Object}{Problem object} 
 \item{model}{
      An object of class \code{\linkS4class{modelorg}}.
    }
    \item{lpdir}{
      Single character string containing the direction of optimization.
      Can be set to \code{"min"} or \code{"max"}.\cr
      Default: \code{"max"}.
    }
    \item{LHS }{constraint matrix, consists of S matrix and ccFBA constraints}
\item{rlb,rub }{bounds for constraints (i.e rows)}
\item{lb,ub }{bounds for variables (i.e columns)}
\item{rtype }{row constraint type}
\item{obj }{index of objective}

    \item{cnames}{
      A character vector giving the variable names. If set to \code{NULL},
      the reaction id's of \code{model} are used.\cr
      Default: \code{NULL}.
    }
    \item{rnames}{
      A character vector giving the constraint names. If set to \code{NULL},
      the metabolite id's of \code{model} are used.\cr
      Default: \code{NULL}.
    }
    \item{pname}{
      A single character string containing a name for the problem object.\cr
      Default: \code{NULL}.
    }
    \item{scaling}{
      Scaling options used to scale the constraint matrix. If set to
      \code{NULL}, no scaling will be performed
      (see \code{\link{scaleProb}}).\cr
      Default: \code{NULL}.
    }
    \item{writeProbToFileName}{
      A single character string containing a file name to which the problem
      object will be written in LP file format.\cr
      Default: \code{NULL}.
    }
    \item{...}{
      Further arguments passed to the initialize method of
      \code{\linkS4class{sysBiolAlg}}. They are \code{solver},
      \code{method} and \code{solverParm}.
    }
%}
%}
}

\author{
  Abdelmoneim Amer Desouki
}



\seealso{
  Constructor function \code{\link{sysBiolAlg}} and
  superclass \code{\linkS4class{sysBiolAlg}}.
}

\keyword{methods}
\keyword{optimize}
