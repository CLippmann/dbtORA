\name{ontologyNr}
\alias{ontologyNr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns for given GOterm the corresponding gene ontology.
}
\description{
Returns for given GOterm (as ID or number) the corresponding gene ontology as number, where 
1 codes for biological process, 
2 for molecular function and 
4 for cellular component. 
If the result is 0, something went wrong.  
}
\usage{
ontologyNr(GOtermNrOrId, Verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{GOtermNrOrId}{
		String vector [1:n] of GOtermIds (like 'GO:0008150') OR Numeric vector [1:n] 
		of GOterm numbers (like 8150).}
	\item{Verbose}{
		Boolean variable. If TRUE, function prints information in GUI window. 
		Default: FALSE}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{OntoNumber}{Numeric vector [1:n] indicating the Gene Ontology that the GOterm belongs to.}
}

\author{
CL
}


