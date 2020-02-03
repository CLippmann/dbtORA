\name{plotGOgraph}
\alias{plotGOgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to draw and save gene ontology DAG.
}
\description{
Function draws and colors the gene ontology DAG of input GO terms depending on input parameters
and saves it as PlotFile.
}
\usage{
plotGOgraph(Adj,GOtermIDs,PlotFile,PlotDirectory=getwd(),
	Significant=rep(1,length(GOtermIDs)),IsHeadline=rep(0,length(GOtermIDs)),
	MarkDetails=TRUE, Overwrite=TRUE, GOtermString=NULL,Remarkable=NULL,Pvalues=NULL,
	NrGenesInTerm=NULL,Expected = NULL, Observed = NULL, Importance = NULL, Up=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Adj}{Adjacency matrix of GO terms where Adj[i,j]==1 iff i is parent of j.}
\item{GOtermIDs}{GO term IDs, e.g. "GO:0008150".}
\item{PlotFile}{Name of the output file. Extension determines the output format, one of png, eps, svg or pdf. Default output format: png.}
\item{PlotDirectory}{Output directory. Default: getwd().}
\item{Significant}{Terms with \code{Significant==1} are drawn in red, others white.}
\item{IsHeadline}{Terms with \code{IsHeadline==1} are marked in yellow.}
\item{MarkDetails}{Set TRUE if details of the dag should be coloured blue. Default: TRUE.}
\item{Overwrite}{Set TRUE if files in \code{PlotDirectory} with same file name as \code{PlotFile} should be replaced by new file. Default: TRUE.}
\item{GOtermString}{GO term strings, e.g "biological_process". If not given (NULL), GO term strings will be generated automatically.}
\item{Remarkable}{Remarkable value of GO terms. See also \code{\link{remarkableness}}.}
\item{Pvalues}{P-values of GO terms. Default: NULL.}
\item{NrGenesInTerm}{The number of annotated genes for each GO term. Default: NULL}
\item{Observed}{Number of observed annotations for each GO term. Default: NULL.}
\item{Expected}{Number of statistically expected annotations for each GO term. Default: NULL.}
\item{Importance}{Importance value of GO terms. Default: NULL.}
\item{Up}{If there are more than the expected number of genes annotated to the GO term, \code{Up} is set to 1 else 0. GO terms where Up==1 are marked in red, others in green. Default: NULL}
}

\value{
PlotFile in PlotDirectory with DAG of input GO terms.
}

\author{
CL
}

\seealso{
For further details about coloring, please see \code{\link{drawORA}}.}


