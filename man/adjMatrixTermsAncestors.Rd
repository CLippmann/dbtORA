\name{adjMatrixTermsAncestors}
\alias{adjMatrixTermsAncestors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to return the adjacency matrix of input GO term numbers in specified ontology.
}
\description{
Function to return the adjacency matrix of input GO term numbers in specified ontology.
}
\usage{
adjMatrixTermsAncestors(GOtermNrInclAncestors, OntologyNr = 1)
}
\arguments{
  \item{GOtermNrInclAncestors}{
Numeric; 
Vector of GO term numbers of GO terms and their 
ancestors up to the root in gene ontology specified 
by \code{OntologyNr}.}
  \item{OntologyNr}{
Numeric; Default: 1
To select the ontology. One of 
1 for biological process, 
2 for molecular function or 
4 for cellular component.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{AdjMatrix}{Numeric;
Adjacency matrix of GOtermNrInclAncestors. 
AdjMatrix[i,j] == 1 iff GO term i is parent of 
GO term j. Named by GOtermIds.}
\item{GOtermNrs}{Numeric;
The GO term numbers corresponding to rows and columns
of AdjMatrix.}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{termsAncestors}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
