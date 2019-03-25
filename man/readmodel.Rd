\name{readmodel}
\alias{readmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
read MOMENT model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
create lp from lists generated from MATLAB MOMENT model.
}
\usage{
readmodel(mat, mets, rxns, rbnds, cbnds, solver = "glpkAPI")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
contain the constraints matrix
}
  \item{mets}{
%%     ~~Describe \code{mets} here~~
list of metabolites
}
  \item{rxns}{
%%     ~~Describe \code{rxns} here~~
list of reactions and their bounds
}
  \item{rbnds}{
%%     ~~Describe \code{rbnds} here~~
bounds of rows of constraint matrix
}
  \item{cbnds}{
%%     ~~Describe \code{cbnds} here~~
bounds of columns of constraint matrix

}
  \item{solver}{
%%     ~~Describe \code{solver} here~~
solver used to solve the lp, can be glpkAPI or cplexAPI
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
return fluxes obtained using the lp.
}
\references{Adadi, R., Volkmer, B., Milo, R., Heinemann, M., & Shlomi, T. (2012). 
Prediction of Microbial Growth Rate versus Biomass Yield by a Metabolic Network 
with Kinetic Parameters, 8(7). doi:10.1371/journal.pcbi.1002575}

\author{Abdelmoneim Amer Desouki}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{cfba_moment}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
