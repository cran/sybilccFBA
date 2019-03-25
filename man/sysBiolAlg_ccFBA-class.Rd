\name{sysBiolAlg_ccFBA-class}

\Rdversion{1.1}
\encoding{utf8}

\docType{class}

\alias{sysBiolAlg_ccFBA-class}
\alias{sysBiolAlg_ccFBA}
\alias{ccFBA}

\title{Class \code{"sysBiolAlg_ccFBA"}}

\description{
  The class \code{sysBiolAlg_ccFBA} holds an object of class
  \code{\linkS4class{optObj}} which is generated to meet the
  requirements of the ccFBA algorithm.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  
  \code{sysBiolAlg(model, algorithm = "ccFBA", ...)}.
  
  Arguments to \code{...} which are passed to method \code{initialize} of class
  \code{sysBiolAlg_ccFBA} are described in the Details section.
}

\section{Slots}{
  \describe{
    \item{\code{problem}:}{Object of class \code{"optObj"}
      containing the problem object.
    }
    \item{\code{algorithm}:}{Object of class \code{"character"}
      containing the name of the algorithm.
    }
    \item{\code{nr}:}{Object of class \code{"integer"}
      containing the number of rows of the problem object.
    }
    \item{\code{nc}:}{Object of class \code{"integer"}
      containing the number of columns of the problem object
    }
    \item{\code{fldind}:}{Object of class \code{"integer"}
      pointers to columns (variables) representing a flux (reaction) in the
      original network. The variable \code{fldind[i]} in the problem object
      represents reaction \code{i} in the original network.
    }
    \item{\code{alg_par}:}{Object of class \code{"list"}
      containing a named list containing algorithm specific parameters.
    }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{sysBiolAlg}"}, directly.
}

\section{Methods}{
  No methods defined with class "sysBiolAlg_ccFBA" in the signature.
}

\details{
  The \code{initialize} method has the following arguments:
  \describe{
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
  }

  The problem object is built to be capable to perform cost constraint flux balance analysis
  (FBA) with a given model, which is basically the solution of a linear
  programming problem
  \deqn{%
      \begin{array}{rll}%
          \max            & \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
                                         \mbox{\boldmath$v$\unboldmath} \\[1ex]
          \mathrm{s.\,t.} & \mbox{\boldmath$Sv$\unboldmath} = 0         \\[1ex]
                          & \alpha_i \leq v_i \leq \beta_i
                            & \quad \forall i \in \{1, \ldots, n\}      \\[1ex]
      \end{array}%
  }{
      max  c^T v
      s.t. Sv = 0
           a_i <= v_i <= b_i  for i = 1, ..., n
  }
  with \eqn{\bold{S}}{S} being the stoichiometric matrix, \eqn{\alpha_i}{a_i}
  and \eqn{\beta_i}{b_i} being the lower and upper bounds for flux (variable)
  \eqn{i} respectively.  The total number of variables of the optimization
  problem is denoted by \eqn{n}.  The solution of the optimization is a flux
  distribution maximizing the objective function
  \eqn{
      \mbox{\boldmath$c$\unboldmath}^{\mathrm{T}}
      \mbox{\boldmath$v$\unboldmath}
  }{c^T v} under the a given environment and the assumption of steady state.
  The optimization can be executed by using \code{\link{optimizeProb}}.
}

\references{
Adadi, R., Volkmer, B., Milo, R., Heinemann, M., & Shlomi, T. (2012). 
Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network 
with Kinetic Parameters, 8(7). doi:10.1371/journal.pcbi.1002575

Gelius-Dietrich, G., Desouki, A. A., Fritzemeier, C. J., & Lercher, M. J. (2013). 
sybil–Efficient constraint-based modelling in R. BMC systems biology, 7(1), 125.

}

\author{
 Abdelmoneim Amer Desouki}


\seealso{
  Constructor function \code{\link{sysBiolAlg}} and
  superclass \code{\linkS4class{sysBiolAlg}}.
}

\keyword{classes}
