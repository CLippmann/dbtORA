\name{checkORAparameters}
\alias{checkORAparameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to check if the parameters passed to \code{dbtORA} are correct.
}
\description{
Function to check if the parameters passed to \code{dbtORA} are correct.
}
\usage{
checkORAparameters(InFileWithExt, InFileDirectory, RefSetFileWithExt, 

RefSetDirectory, OutFile, OutFileDirectory, Correction, PvalueThreshold,

MinNrOfGenes, OnlyManuCur, drawDAG, MarkDetails, MarkHeadlines, PlotExt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InFileWithExt}{
String; 
Filename with extension where NCBIs are keys (for *.names and *.lrn files)
or the only column (for *.txt files).}
  \item{InFileDirectory}{
String;
Directory where InFileWithExt can be found.}
  \item{RefSetFileWithExt}{
String;
Filename with extension where NCBIs are keys (for *.names and *.lrn files)
or the only column (for *.txt files). NCBIs will be used as reference set.
}
  \item{RefSetDirectory}{
String;
Directory where RefSetFileWithExt with reference NCBIs can be found.
}
  \item{OutFile}{
String;
Filename of the output file(s). Will be complemented by the parameters of the ORA.
}
  \item{OutFileDirectory}{
String;
Directory where results of ORA and DAGs will be saved.
}
  \item{Correction}{
String;
Type of correction for mulitple testing of the p-values.
'BON' for Bonferroni, 
'FDR' for False Discovery Rate,
'RAW' if no correction should be done.}
  \item{PvalueThreshold}{
Numeric;
P-value threshold. GO-Terms with p-values greater than PvalueThreshold will be ignored.
}
  \item{MinNrOfGenes}{
Numeric;
Minimum number of genes annotated to one Term that is accepted. GO-Terms with more than
MinNrOfGenes annotated genes will be considered in calculation.
}
  \item{OnlyManuCur}{
Boolean;
Set TRUE if only manually curated gene annotations should be considered.
}
  \item{drawDAG}{
Boolean;
Set TRUE if directed acyclic graphs (DAGs) should be drawn.
If drawDAG is set to FALSE, the parameters MarkDetails, MarkHeadlines and CurtDAG will be
ignored.
}
  \item{MarkDetails}{
Boolean;
Set TRUE if details of the DAG should be marked in blue colour.
}
  \item{MarkHeadlines}{
Boolean;
Set TRUE if headlines of the DAG should be marked in yellow colour.
}
  \item{PlotExt}{
String;
Extension of the plotfile showing the DAG. One of 'pdf' , 'eps' or 'png'.
}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dbtORA}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
