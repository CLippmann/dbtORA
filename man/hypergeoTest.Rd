\name{hypergeoTest}
\alias{hypergeoTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do a one-sided hypergeometric test.
}
\description{
Function to do a one-sided hypergeometric test, i.e. calculate the probability to draw more or less 
(expectation value smaller than observed number of successes respectively expectation value greater 
than observed number of successes) than a certain number of successes (ObservedNrOfAnnsInTerm) in a 
fixed number of draws (NrOfGenesInSample), without replacement, from a finite population of fixed 
size (NrOfGenesInUniverse) that contains a known number of successes (NrOfAnnotationsInTerm), wherein 
each draw is either a success or a failure.
}
\usage{
hypergeoTest(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, 

NrOfGenesInUniverse, LogPvalues = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ObservedNrOfAnnsInTerm}{
Numeric;
Vector of observed numbers of input genes annotated to one GO term. }
  \item{NrOfAnnotationsInTerm}{
Numeric;
Vector of numbers of all genes annotated to one GO term.}
  \item{NrOfGenesInSample}{
Numeric;
The number of input genes (genes of interest in sample) annotated to 
at least one GO term.}
  \item{NrOfGenesInUniverse}{
Numeric;
The number of genes in universe, i.e. all genes annotated to at 
least one GO term.}
  \item{LogPvalues}{
Boolean; Default: TRUE
Set TRUE if -log(p-values) should be calculated. Set FALSE if non-transformed p-values should be returned.
}
}
\details{
Hypergeometric test is done one-sided depending on \code{ExpectedNrOfAnnsInTerm}:
If the expected number of genes annotated to one GO term is less than 
ObservedNrOfAnnsInTerm, the log-p-value will be log(P(X>=ObservedNrOfAnnsInTerm)) 
where X is the hypergeometric distributed random variable. 
If the expected number of genes annotated to one GO term is greater than 
ObservedNrOfAnnsInTerm, thelog-p-value will be log(P(X<ObservedNrOfAnnsInTerm)) 
where X is the hypergeometric distributed random variable.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{LogPvalues}{Numeric;
Vector of log-p-values of one-sided hypergeometric test.}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{phyper}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
