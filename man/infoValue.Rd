\name{infoValue}
\alias{infoValue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function calculates the partial Shannon information of gene sets in GO terms.
}
\description{
Function calculates the partial Shannon information of gene sets in GO terms
explaining how informative a certain term in the context of all terms is.
}
\usage{
infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse = max(NrOfAnnotationsInTerm))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NrOfAnnotationsInTerm}{
Numeric; 
Vector of numbers of genes annotated to corresponding GO terms.}
  \item{NrOfGenesInUniverse}{
Numeric; Default: max(NrOfAnnotationsInTerm)
Number of genes in universe. (If not restricted to reference set, NrOfGenesInUniverse
is the same as the number of genes (directly + indirectly) annotated to the root.)}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{InfoValue}{Numeric;
A value for each term that describes how informative that term is.}
\item{InfoValueP}{Numeric;
The ratio of NrOfAnnotationsInTerm/NrOfGenesInUniverse for each term which is
the empirical probability of occurrence.}
}

\author{
CL
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{certainty}}, 
\code{\link{importance}}, 
\code{\link{remarkableness}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
