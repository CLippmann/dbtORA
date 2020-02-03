\name{NCBI2GeneName}
\alias{NCBI2GeneName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to get GeneSymbol and GeneName for given NCBI numbers.
}
\description{Function to get GeneSymbol and GeneName for given NCBI numbers.
From ReadNAMES('AllAnnNCBIsPlusGeneName.names', system.file('extdata',package='ORA')).}
\usage{
NCBI2GeneName(NCBI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NCBI}{Numeric, vector of NCBI numbers (e.g. 1).}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{GeneSymbol}{String, vector of same length as NCBI containing abbreviation
of GeneName (e.g. "A1BG").}
\item{GeneName}{String, vector of same length as NCBI containing detailed
description of the gene (e.g. "alpha-1-B glycoprotein").}
}

\author{
CL
}



