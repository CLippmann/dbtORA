\name{drawORA}
\alias{drawORA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To draw the gene ontology DAG containing ORA results' information.
}
\description{
To draw the gene ontology DAG containing ORA results' information.
Prepares data for \code{plotGOgraph} which does the actual plotting.
}
\usage{
drawORA(ORAresults, PlotFileWithExt, PlotDirectory, MarkDetails = TRUE, 

MarkHeadlines = TRUE, Overwrite = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ORAresults}{
		List of 4:
		\describe{
			\item{LRNresults}{List of 16:
				Information needed to generate the lrn file containing all the calculated values for
				the GO terms found to be significant for input genes.}
			\item{NAMESresults}{List of 3:
				Information needed to generate the names file containing information about GO terms.}
			\item{Genes2GOtermsSparseMatrix}{Numeric; Sparse matrix explaining the connection of genes and GO terms. 
				Genes2GOtermsSparseMatrix[i,3]==1	iff gene in ith row is annotated to 
				GO term in ith row.}
			\item{GO2GOAdjMatrices}{List of 4:
				Adjacency matrices for each ontology and combined sparse matrix.}
		}
	}
  \item{PlotFileWithExt}{
		String;
		Name of the file that should be drawn with extension.
		Extension can be one of 'png', 'eps' or 'pdf'. Default is 'png'.}
  \item{PlotDirectory}{
		String; Default: current directory
		Directory where PlotFileWithExt should be saved.}
  \item{MarkDetails}{
		Boolean; Default: TRUE
		Set TRUE if details of DAG should be marked in blue colour.}
  \item{MarkHeadlines}{
		Boolean; Default: TRUE
		Set TRUE if headlines should be marked in yellow colour in DAG.}
  \item{Overwrite}{
		Boolean; Default: TRUE
		Set TRUE if existing files with the same name should be overwritten.}
}

\details{	Coloring of the nodes and its meaning:
					Significantly overrepresented nodes are marked in red; 
					Significantly underrepresented nodes are marked in green;
					Terms that are important for DAG structure but do not have a
					significant p-value are left white. 
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

\author{
CL
}

\seealso{
For further details about ORAresults, please see \code{ORA}.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
