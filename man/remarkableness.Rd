\name{remarkableness}
\alias{remarkableness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the remarkableness of a GO term.
}
\description{
Function to calculate the remarkableness of a GO term.
}
\usage{
remarkableness(Certainty, InfoValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Certainty}{
Numeric;
(Vector of) certainty value(s) of GO term.
See certainty(Pvalues).}
  \item{InfoValue}{
Numeric;
(Vector of) information value(s) of GO term.
See infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse).}
}

\value{
\item{Remakable}{Numeric;
Vector containing remarkableness values = certainty * information value of the GO term(s).}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{certainty}}, 
\code{\link{infoValue}}, 
\code{\link{importance}}.
}


