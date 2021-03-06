\name{addGlcTrns}
\alias{addGlcTrns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
add Glucose Transport constraint }
\description{
add glucose transport constraint to the problem. Put an upperbound on glucose consumption.
}
\usage{
addGlcTrns(prob, mod2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
lp problem
  }
  \item{mod2}{
 An object of class \code{\link{modelorg}} with only irreversible reactions.
        It can be sent  to save time of recalculating it with each call.
  }
}
\author{Abdelmoneim Amer Desouki}

\seealso{
\code{\link{modelorg}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (prob, mod2) 
{
    Si = 0.2
    Hxt1 = 41 * Si/(Si + 107)
    Hxt2 = 16.1 * Si/(Si + 2.9)
    Hxt3 = 18.5 * Si/(Si + 29)
    Hxt4 = 12 * Si/(Si + 6.2)
    Hxt5 = 14 * Si/(Si + 10)
    Hxt6 = 11.4 * Si/(Si + 1.5)
    Hxt7 = 11.7 * Si/(Si + 1.3)
    Gal2 = 17.5 * Si/(Si + 1.5)
    colid = getNumCols(lp = problem(prob)) + 1
    trnsCol = NULL
    rowind = getNumRows(lp = problem(prob)) + 1
    glcRxn = which(react_id(mod2) == "R_GLCt1")
    addRowsToProb(lp = problem(prob), i = rowind, type = "U", 
        lb = 0, ub = 0, cind = list(c(trnsCol[1, "Col"], trnsCol[2, 
            "Col"], trnsCol[3, "Col"], trnsCol[4, "Col"], trnsCol[5, 
            "Col"], trnsCol[6, "Col"], trnsCol[7, "Col"], glcRxn)), 
        nzval = list(c(-Hxt1, -Hxt2, -Hxt3, -Hxt4, -Hxt5, -Hxt6, 
            -Hxt7, 1)), rnames = "glcTrns")
    return(prob)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
