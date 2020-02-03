\name{termsAncestors}
\alias{termsAncestors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function returns vector of all ancestors in ontology for given GO term numbers.
}
\description{
Function returns vector of all ancestors in ontology for given GO term numbers.
}
\usage{
termsAncestors(GOtermNr, OntologyNr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOtermNr}{
Numeric; 
Vector of GO term numbers.
}
  \item{OntologyNr}{
Numeric;
To select the ontology. One of: 
1 for biological process, 
2 for molecular function or 
4 for cellular component
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Ancestors}{Numeric; Unique and by number sorted GO terms that are ancestors of input GO terms.}
\item{TermsWithoutAncestors}{Numeric; Vector of those GO terms for that no ancestors could be found.}
}

\author{
CL
}


