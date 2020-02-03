\name{ORAfilename}
\alias{ORAfilename}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to complete the OutFile name passed to dbtORA with ORA parameters.
}
\description{
Function to complete the OutFile name passed to dbtORA with ORA parameters.
}
\usage{
ORAfilename(OutFile, NrOfValidInputGenes, Correction, PvalueThreshold, MinNrOfGenes, 
OnlyManuCur, WithRefSet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OutFile}{
String;
Filename of the output file(s). Will be complemented by the parameters of the ORA.}
	\item{NrOfValidInputGenes}{
String;
Number of valid input genes = #(input genes) - #(duplicated and non-annotated genes)
	}
  \item{Correction}{
String;
Type of correction for mulitple testing of the p-values.
'BON' for Bonferroni, 
'FDR' for False Discovery Rate,
'RAW' if no correction should be done.}
  \item{PvalueThreshold}{
Numeric;
P-value threshold. GO-Terms with p-values greater than PvalueThreshold will be ignored.}
  \item{MinNrOfGenes}{
Numeric;
Minimum number of genes annotated to one Term that is accepted. Only GO terms with more than
MinNrOfGenes genes will be considered in calculation.}
  \item{OnlyManuCur}{
Boolean;
Set TRUE if only manually curated gene annotations should be considered.}
	\item{WithRefSet}{
	Boolean; Default: FALSE. Set TRUE if a reference set of genes is used.
	}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{OutFilePlusParams}{String;
The complemented OutFile.}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dbtORA}}
}


