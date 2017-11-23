adjMatrixTermsAncestors <- function(GOtermNrInclAncestors, OntologyNr = 1){
# Function to return the adjacency matrix of input GOterm numbers in specified
# ontology.

# AdjMatrixGOterms2GOterms <- adjMatrixTermsAncestors(GOtermNrInclAncestors, OntologyNr)

# INPUT:
# GOtermNrInclAncestors		Numeric; 
#													Vector of GO term numbers of GO terms and their 
#													ancestors up to the root in gene ontology specified 
#													by OntologyNr. 
#
# OPTIONAL:
# OntologyNr							Numeric; Default: 1
#													To select the ontology. One of 
#													1 for biological process, 
#													2 for molecular function or 
#													4 for cellular component.
#
# OUTPUT:
# AdjMatrix								Numeric;
#													Adjacency matrix of GOtermNrInclAncestors. 
#													AdjMatrix[i,j] == 1 iff GO term i is parent of 
#													GO term j. Named by GOtermIds.
# GOtermNrs								Numeric;
#													The GO term numbers corresponding to rows and columns
#													of AdjMatrix.

# EXAMPLE:
# GOtermeMitAncestors <- union(termsAncestors(c(1,2,3,11,22), 1)$Ancestors, c(1,2,3,11,22))
# OntologyNr <- 1
# AdjMatrixGOterms2GOterms <- adjMatrixTermsAncestors(GOtermeMitAncestors, OntologyNr)

# AUTHOR:
# CL, 30.03.2016

# USES:
# packages: 'Matrix'
# functions:
 # [1] "function"         "require"          "if"               "length"          
 # [5] "stop"             "all"              "is.numeric"       "switch"          
 # [9] "as.character"     "ReadSparseMatrix" "GOdataDi"         "rownames"        
# [13] "match"            "termNr"           "any"              "is.na"           
# [17] "which"            "as.logical"       "rowSums"          "sapply"          
# [21] "union"            "matrix"           "sort"             "colnames"        
# [25] "cbind"            "as.numeric"       "return"           "list"            


#require('Matrix')
#requireNamespace(package ='Matrix', quietly = TRUE)

# Dateneingabe ueberpruefen:
if(length(GOtermNrInclAncestors) <= 0){stop('adjMatrixTermsAncestors: A non-empty vector of GO term numbers including the terms ancestors is required as input! Function stops.')}
if(!all(is.numeric(GOtermNrInclAncestors))){stop('adjMatrixTermsAncestors: GO term numbers are not numeric! Function stops.')}
if(!is.numeric(OntologyNr)){stop('adjMatrixTermsAncestors: OntologyNr has to be numeric: 1 (=BP), 2 (=MF) or 4 (=CC)! Function stops.')}


# kinderadjmatrix einlesen je nach OntologyNr:
switch(as.character(OntologyNr),
	'1' =	{	# biological process = BP
					kinderadjsparse <- ReadSparseMatrix('AdjBPsparseMatrix', GOdataDi('09Originale'))
					GOtermNames <- kinderadjsparse$DimNames$rownames
				},
	'2' = { # molecular function = MF
					kinderadjsparse <- ReadSparseMatrix('AdjMFsparseMatrix', GOdataDi('09Originale'))
					GOtermNames <- kinderadjsparse$DimNames$rownames
				},
	'4' = { # cellular component = CC
					kinderadjsparse <- ReadSparseMatrix('AdjCCsparseMatrix', GOdataDi('09Originale'))
					GOtermNames <- kinderadjsparse$DimNames$rownames
				},
	stop('termsAncestors: OntologyNr has to be 1 (=BP), 2 (=MF) or 4 (=CC)! Function stops.')
)# end switch(OntologyNr)


# terme matchen
termeInd <- match(GOtermNrInclAncestors, termNr(GOtermNames))
if(any(is.na(termeInd))){stop("adjMatrixTermsAncestors: Some GO term numbers couldn't be found in specified gene ontology. A vector of GO term numbers including the terms' ancestors that are all in the same, specified gene ontology, is required as input! Function stops.")}
SpaltenIndIND <- which(as.logical(Matrix::rowSums(sapply(termeInd,'==',kinderadjsparse$SpaltenInd))))

# nur die Infos zu GOtermNrInclAncestors selektieren:
ZeilenInd <- kinderadjsparse$ZeilenInd[SpaltenIndIND]
SpaltenInd <- kinderadjsparse$SpaltenInd[SpaltenIndIND]
Inhalt <- kinderadjsparse$Inhalt[SpaltenIndIND]

# jetzt ne kleinere ADJ-Matrix (nicht sparse!) daraus machen: 
VorkommendeIndizes <- union(ZeilenInd, SpaltenInd)
Dim <- length(VorkommendeIndizes) # Dimension
if(Dim==0){ # Wenn die Dimension Null ist, leere Adj zurueckgeben
	if(length(termeInd)==1){ # Spezialfall: Es ist nur die Wurzel in der Adj also length(termeInd)==1 dann ne Adj mit nur der Wurzel zurueckgeben.
		AdjMatrix <- matrix(0, ncol=1, nrow=1)
		rownames(AdjMatrix) <- GOtermNames[termeInd]
		colnames(AdjMatrix) <- GOtermNames[termeInd]
		return(list(AdjMatrix = AdjMatrix, GOtermNrs = termNr(GOtermNames)[termeInd]))
	}
	return(list(AdjMatrix = matrix(0, ncol=Dim, nrow=Dim), GOtermNrs = NULL))
}
AdjMatrix <- matrix(0, ncol=Dim, nrow=Dim)
rownames(AdjMatrix) <- sort(VorkommendeIndizes)
colnames(AdjMatrix) <- sort(VorkommendeIndizes)
AdjMatrix[cbind(as.character(ZeilenInd), as.character(SpaltenInd))] <- 1


# GOterme finden, die zu den row- and colnames der AdjMatrix gehoeren:
GOtermIds <- GOtermNames[as.numeric(rownames(AdjMatrix))]
rownames(AdjMatrix) <- GOtermIds
colnames(AdjMatrix) <- GOtermIds

return(list(AdjMatrix = AdjMatrix, GOtermNrs = termNr(GOtermIds)))
}# end function adjMatrixTermsAncestors
