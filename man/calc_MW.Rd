\name{calc_MW}
\alias{calc_MW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate molecular weights
}
\description{
Calculate Molecular weights of different proteins using the genome .faa file.
}
\usage{
calc_MW(aa_fname = "aa.txt", ptt_fname = "test2.ptt", faa_fname = "NC_000913.faa",
 nchrm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aa_fname}{
file name of file containing lisT of amino acid names
}
  \item{ptt_fname}{
file name of file containing gene names with gene code
}
  \item{faa_fname}{
file name of file containing gene code and sequence of amino acids
}
 \item{nchrm}{
the number of chromosomes in the genome
}
}

\value{
generate a file containing gene name , length, and molecular weight

}

\author{Abdelmoneim Amer Desouki}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\dontrun{
aa_fname <- system.file("extdata", "aa.txt", package="sybilccFBA")
ptt_fname <- system.file("extdata", "test2.ptt", package="sybilccFBA")
faa_fname <- system.file("extdata", "NC_000913.faa", package="sybilccFBA")

geneCnt <- calc_MW(aa_fname,ptt_fname,faa_fname)
write.csv(file="geneCnt.csv",geneCnt) 

## The function is currently defined as
"calc_MW"
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Molecular weights }
