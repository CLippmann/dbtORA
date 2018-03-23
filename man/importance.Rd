\name{importance}
\alias{importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the importance for given information values and certainty values for GO terms.
}
\description{
Function to calculate the importance for given information values and certainty values for GO terms.
}
\usage{
importance(Certainty, InfoValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Certainty}{
Numeric; Values between 0 and 100.
Vector of certainty values for GO terms. See certainty(Pvalues).}
  \item{InfoValue}{
Numeric; Values between 0 and 100.
Vector of (partial shannon) information values for GO terms. 
See infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse).}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Importance}{Numeric; Values between 0 and 100.
Vector of importance value based on InfoValue and Certainty, i.e. minimum of InfoValue and Certainty.}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{certainty}}, 
\code{\link{infoValue}}, 
\code{\link{remarkableness}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
