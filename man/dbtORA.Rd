\name{dbtORA}
\alias{dbtORA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to do an ORA including the drawing of DAGs.
}
\description{
Function to do an Overrepresentation Analysis (ORA) including the drawing of the directed acyclic graphs (DAGs)
of the resulting GO terms.
}
\usage{
dbtORA(InFileWithExt, PvalueThreshold = 0.05, Correction = "BON", OnlyManuCur = TRUE, 

MinNrOfGenes = 2, InFileDirectory = getwd(), OutFile = InFileWithExt, 

OutFileDirectory = InFileDirectory, RefSetFileWithExt = NULL, 

RefSetDirectory = InFileDirectory, drawDAG = TRUE, MarkDetails = TRUE, 

MarkHeadlines = TRUE, PlotExt = "png")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InFileWithExt}{
String; 
Filename with extension where NCBIs are keys (for *.names and *.lrn files)
or the only column (for *.txt files).}
  \item{PvalueThreshold}{
Numeric; Default: 0.05
P-value threshold. GO terms with p-values greater than PvalueThreshold will be ignored.}
  \item{Correction}{
String; Default: 'BON'
Type of correction for mulitple testing of the p-values.
'BON' for Bonferroni, 
'FDR' for False Discovery Rate,
'RAW' if no correction should be done.}
  \item{OnlyManuCur}{
Boolean; Default: TRUE
Set TRUE if only manually curated gene annotations should be considered.}
  \item{MinNrOfGenes}{
Numeric; Default: 2
Minimum number of genes annotated to one Term that is accepted. Only GO terms with more than
MinNrOfGenes annotated genes will be considered in calculation.}
	\item{InFileDirectory}{
String; Default: current directory getwd().
Directory where InFileWithExt can be found.}
  \item{OutFile}{
String; Default: InFileWithExt (extension will be adjusted)
Filename of the output file(s). Will be complemented by the parameters of the ORA.}
  \item{OutFileDirectory}{
String; Default: InFileDirectory
Directory where results of ORA and DAGs will be saved.}
  \item{RefSetFileWithExt}{
String; Default: NULL
Filename with extension where NCBIs are keys (for *.names and *.lrn files)
or the only column (for *.txt files). NCBIs will be used as reference set.}
  \item{RefSetDirectory}{
String; Default: InFileDirectory
Directory where RefSetFileWithExt with reference NCBIs can be found.}
  \item{drawDAG}{
Boolean; Default: TRUE
Set TRUE if DAGs should be drawn.
If drawDAG is set to FALSE, the parameters MarkDetails, MarkHeadlines and CurtDAG will be
ignored.}
  \item{MarkDetails}{
Boolean; Default: TRUE
Set TRUE if details of the DAG should be marked in blue colour.}
  \item{MarkHeadlines}{
Boolean; Default: TRUE
Set TRUE if headlines of the DAG should be marked in yellow colour.}
  \item{PlotExt}{
String; Default: 'png' 
Extension of the plotfile showing the DAG. One of 'pdf', 'eps' or 'png'.
}
}
\details{
	Wrapper function to execute mainly \code{ORA} and \code{drawORA}.
	Coloring of the nodes and its meaning:
	Red - Significantly overrepresented nodes; 
	Green - Significantly underrepresented nodes;
	White - Terms that are important for DAG structure but do not have a
	significant p-value. 
	If MarkHeadlines = TRUE, the significant(!) nodes with highest remarkable value
	in each path from a detail to the root get a yellow filling. The 
	margin indicates over- or underrepesentation by its red or green color.
	If MarkDetails = TRUE, the details of the DAG will be colored in blue. 
	The margin again indicates over- or underrepesentation by its red or green color.
	If MarkHeadlines and MarkDetails are TRUE, there might be nodes that are
	both headlines and details. In this case the nodes have a maring acording to
	over- or underrepesentation in red or green and are filled in yellow like all
	headlines. Additionally the writing is blue to indicate that this node is a detail.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Various files containing all the information received by ORA, Genes to GO terms matrix, adjacency matrix of the 
GO terms and plots of the corresponding DAGs.
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ORA}},
\code{\link{drawORA}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
