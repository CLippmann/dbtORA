\name{WriteORAresults}
\alias{WriteORAresults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to write lrn file of GO terms and the computed values, names file of GO terms and their
descriptions and sparse matrix of genes and GO terms where the results of the overrepresentation
analysis (ORA) are given.
}
\description{
Function to write lrn file of GO terms and the computed values, names file of GO terms and their
descriptions and sparse matrix of genes and GO terms where the results of the overrepresentation
analysis (ORA) are given.
}
\usage{
WriteORAresults(FileNameWithoutExt, ORAresults, OutDirectory = getwd(), 

InFileWithExt = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FileNameWithoutExt}{
			String.
			Name of the output file without extention.}
  \item{ORAresults}{
			List of 4. For further details see function ORA.
			\describe{
					\item{LRNresults}{List of 19 containing results relevant for lrn:
							LRNresults = list(GOtermNr, OntologyNr, NrOfGenesInUniverse, 
							NrOfGenesInSample, NrOfAnnotationsInTerm, Up, 
							ExpNrOfAnnsInTerm, ObservedNrOfAnnsInTerm, 
							RelDiff, Pvalue, LogPvalue, Certainty, InfoValue, 
							Remarkable, Importance, InfoContent, InfoContentORA, IsHeadline, IsDetail)
							Names of this list will be used as Header in FileNameWithoutExt.lrn.}
					\item{NAMESresults}{List of 3 containing results relevant for names:
							NAMESresults = list(GOtermNr, GOtermDescription, GOtermId)
							Names of this list will be used as Header in FileNameWithoutExt.names.}	
					\item{Genes2GOtermsMatrix}{Matrix describing the annotations of genes to GO terms:
							Dimnames will be used as Key (rownames) respectively Header (colnames).}
					\item{GO2GOAdjMatrices}{List of 4 containing the adjacency matrices for each ontology and the 
							combined sparse matrix representing the DAGs of significant GOterms:
							GO2GOAdjMatrices = list(GO2GOSparseAdjMatrix, AdjMatrixGO2GOBP, 
												AdjMatrixGO2GOMF, AdjMatrixGO2GOCC)}
			}
	}	
	\item{OutDirectory}{
			String. Default: current directory.
			Directory where the files should be saved.
	}
	\item{InFileWithExt}{
			String. Default: ''.
			FileName of original input file.
	}
}

\value{
Files saved in OutDirectory containing all the information recieved by ORA, Genes to GO terms matrix
and adjacency matrix of the GO terms.
}

\author{
CL
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
