\name{ORA}
\alias{ORA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate a overrepresentation analysis.
}
\description{
Function to calculate a overrepresentation analysis based on gene ontology
using a one-sided hypergeometric test for given genes.
}
\usage{
ORA(NCBIs, Correction = "BON", PvalueThreshold = 0.05, MinNrOfGenes = 2, 

OnlyManuCur = FALSE, RefSet = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NCBIs}{
Numeric;
Vector of NCBI numbers of genes in sample (gene set of interest).}
  \item{Correction}{
String; Default: 'BON'
Type of correction for mulitple testing of the p-values.
'BON' for Bonferroni, 
'FDR' for False Discovery Rate,
'RAW' if no correction should be done.}
  \item{PvalueThreshold}{
Numeric; Default: 0.05
P-value threshold. GO terms with p-values greater than PvalueThreshold will be ignored.}
  \item{MinNrOfGenes}{
Numeric; Default: 2
Minimum number of genes annotated to one term that is accepted. Only GO terms with more than
MinNrOfGenes genes will be considered in calculation.}
  \item{OnlyManuCur}{
Boolean; Default: TRUE
Set TRUE if only manually curated gene annotations should be considered.}
  \item{RefSet}{
Numeric; Default: NULL
Vector of NCBI numbers of genes that form the reference set (universe). If not given i.e.
NULL or missing, all known genes are taken as universe.}
}

\value{
\item{ORAresults}{
	List of 4:
	\describe{
		\item{LRNresults}{List of 16:
			Information needed to generate the lrn file containing all the calculated values for
			the GO terms found to be significant for input genes.
			\describe{
				\item{LRNresults$GOtermNr}{Numeric; GO term numbers found to be significant for input genes.} 
				\item{LRNresults$OntologyNr}{Numeric; Number of ontology (1=BP, 2=MF, 4=CC).} 
				\item{LRNresults$NrOfGenesInUniverse}{Numeric; Number of genes in universe used for p-value computation.} 
				\item{LRNresults$NrOfGenesInSample}{Numeric; Number of input genes used for p-value computation.} 
				\item{LRNresults$NrOfAnnotationsInTerm}{Numeric; Number of annotations associated to GO term.} 
				\item{LRNresults$Up}{Numeric; 1 if GO term is up regulated (expected < observed), 0 if down.} 
				\item{LRNresults$ExpNrOfAnnsInTerm}{Numeric; Statistically expected number of genes annotated to GO term.} 
				\item{LRNresults$ObservedNrOfAnnsInTerm}{Numeric; Empirically observed number of genes annotated to GO term.}
				\item{LRNresults$RelDiff}{Numeric; Relative difference of expected and observed in percent.}
				\item{LRNresults$Pvalue}{Numeric; Pvalues for each GO term recieved by statistical test.}
				\item{LRNresults$LogPvalue}{Numeric; log(Pvalues).} 
				\item{LRNresults$Certainty}{Numeric; Certainty value. See \code{certainty}.} 
				\item{LRNresults$InfoValue}{Numeric; Value describing partial Shannon information. See \code{infoValue}.}
				\item{LRNresults$Remarkable}{Numeric; Product of Certainty and InfoValue.}
				\item{LRNresults$Importance}{Numeric; Minimum of Certainty and InfoValue.}
				\item{LRNresults$InfoContent}{Numeric; InformationContent from GOTermInfosBP/MF/CC.lrn depending on OnlyManuCur.}
				\item{LRNresults$InfoContentORA}{Numeric; -log2(ObservedNrOfAnnsInTerm/NrOfGenesInSample)}
				\item{LRNresults$IsHeadline}{Boolean; 1 if GO term is headline, 0 if not.}
				\item{LRNresults$IsDetail}{Boolean; 1 if GO term is detail, 0 if not.}
			}
		}
		\item{NAMESresults}{List of 3:
			Information needed to generate the names file containing information about GO terms.
			\describe{
				\item{NAMESresults$GOtermNr}{Numeric; GO term numbers found to be significant for input genes.}
				\item{NAMESresults$GOtermDescription}{String; Description of GO terms = termDescription(GOtermId).} 
				\item{NAMESresults$GOtermId}{String; GO term Id = termId(GOtermNr).}
			}
		}
		\item{Genes2GOtermsMatrix}{Numeric; matrix explaining the connection of genes and GO terms. 
			Genes2GOtermsMatrix[i,j]==1 iff gene in ith row is annotated to 
			GO term in jth row.}
		\item{GO2GOAdjMatrices}{List of 4:
			Adjacency matrices for each ontology and combined sparse matrix.
			\describe{
				\item{GO2GOAdjMatrices$GO2GOSparseAdjMatrix}{Numeric; Sparse adjacency matrix describing the complete directed 
					acyclic graph (DAG) of the significant GOterms up to the root, 
					i.e. the edges between GOterms and their parents. Dimnames of 
					columns are specifying the ontology - BP/MF/CC.}
				\item{GO2GOAdjMatrices$AdjMatrixGO2GOBP}{Numeric; (non-sparse) Adjacency matrix of BP-DAG. 
					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.}
				\item{GO2GOAdjMatrices$AdjMatrixGO2GOMF}{Numeric; (non-sparse) Adjacency matrix of MF-DAG. 
					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.}
				\item{GO2GOAdjMatrices$AdjMatrixGO2GOCC}{Numeric; (non-sparse) Adjacency matrix of CC-DAG. 
					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.}
			}
		}
	}
}
}

\author{
CL
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
