\name{certainty}
\alias{certainty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate the certainty for a GO term.
}
\description{
Function to calculate the certainty for a GO term, i.e. 
the probability that there is no term with a smaller p-value in the given 
GO subtree.}
\usage{
certainty(Pvalues)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pvalues}{
Numeric;
P-values of the GO terms.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Certainty}{Numeric;
Empiric probability that there is no term with a smaller p-value in the given GO subtree.
}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{infoValue}}, 
\code{\link{importance}}, 
\code{\link{remarkableness}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
