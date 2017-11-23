\name{drawFunctionalAreas}
\alias{drawFunctionalAreas}
\title{
Function to draw DAGs for each functional area where ORA results are given.
}
\description{
Function to draw DAGs for each functional area with functional area as root 
where ORA results and functional area nodes are given.
}
\usage{
drawFunctionalAreas(FAfileName, FAdir, ORAfileName, AdjFileName, ORAdir = FAdir, 

PlotDir = FAdir, Ontology = 'BP', DrawRoot2FAs = FALSE, MarkHeadlines=TRUE, 

MarkDetails=TRUE, PlotFileExt='png')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FAfileName}{String. Filename of the *.names-File that contains the Functional Areas. If not given, filenames (\code{FAfileName}, \code{ORAfileName}, \code{OntoADJfileName}) and directories (\code{FAdir}, \code{ORAdir}) will be requested interactively.}
  \item{FAdir}{String. Directory where \code{FAfileName} can be found.}
  \item{ORAfileName}{String. Filename of the *.lrn-File containing the ORA results.}
  \item{AdjFileName}{String. Filename of the *.lrn-File containing the complete adjacency matrix of GOterms resultig from ORA.}
  \item{ORAdir}{String. Directory where ORAfileName and AdjFileName can be found. Default: FAdir.}
  \item{PlotDir}{String. Directory where the plots (DAGs) should be saved. Default: FAdir.}
  \item{Ontology}{String. Ontology (one of: 'BP', 'CC', 'MF') for which the DAGs of FAs should be drawn. Default: 'BP'.}
  \item{DrawRoot2FAs}{Boolean. If DAG from Ontology root to FAs (in \code{Ontology}) should be drawn, set TRUE. Default: FALSE.}
  \item{MarkHeadlines}{Boolean. If Headlines should be marked in yellow, set TRUE. Default: TRUE.}
  \item{MarkDetails}{Boolean. If Details should be marked in blue, set TRUE. Default: TRUE.}
  \item{PlotFileExt}{String. Extension of the PlotFile. One of "png" or "eps". Default: 'png'.}
}
\details{
Functional Areas are nodes of a DAG known to be important. %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Depending on the number of given functional areas files containing the DAGs from functional areas down to details are created. They are named by the GO term strings of the functional areas and a subsequent number.
}

\author{
CL
}

\seealso{
\code{\link{dbtORA}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
