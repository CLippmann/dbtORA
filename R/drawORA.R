drawORA <- function(ORAresults, PlotFileWithExt, PlotDirectory, MarkDetails = TRUE, MarkHeadlines = TRUE, Overwrite=TRUE){
# To draw the gene ontology DAG containing ORA results' information.

# drawORA(ORAresults, PlotFileWithExt, PlotDirectory, MarkDetails = TRUE, MarkHeadlines = TRUE, Overwrite=TRUE)

# INPUT:
# ORAresults	list of 4:
#	LRNresults					List of 16:
#											Information needed to generate the lrn file containing all the calculated values for
#											the GO terms found to be significant for input genes.
	#	LRNresults$GOtermNr								Numeric; GO term numbers found to be significant for input genes. 
	# LRNresults$OntologyNr							Numeric; Number of ontology (1=BP, 2=MF, 4=CC). 
	# LRNresults$NrOfGenesInUniverse		Numeric; Number of genes in universe used for p-value computation. 
	# LRNresults$NrOfGenesInSample			Numeric; Number of input genes used for p-value computation. 
	# LRNresults$NrOfAnnotationsInTerm	Numeric; Number of annotations associated to GO term. 
	# LRNresults$Up											Numeric; 1 if GO term is up regulated (expected < observed), 0 if down. 
	# LRNresults$ExpNrOfAnnsInTerm			Numeric; Statistically expected number of genes annotated to GO term. 
	# LRNresults$ObservedNrOfAnnsInTerm	Numeric; Empirically observed number of genes annotated to GO term.
	# LRNresults$RelDiff								Numeric; Relative difference of expected and observed in percent.
	# LRNresults$Pvalue									Numeric; Pvalues for each GO term recieved by statistical test.
	# LRNresults$LogPvalue							Numeric; log(Pvalues). 
	# LRNresults$Certainty							Numeric; Certainty value. See function certainty. 
	# LRNresults$InfoValue							Numeric; Value describing partial Shannon information. See infoValue.
	# LRNresults$Remarkable							Numeric; Product of Certainty and InfoValue. 
	# LRNresults$Importance							Numeric; Minimum of Certainty and InfoValue. 
	# LRNresults$IsHeadline							Boolean; 1 if GO term is headline, 0 if not.
	# LRNresults$IsDetail								Boolean; 1 if GO term is detail, 0 if not.
#
# NAMESresults				List of 3:
#											Information needed to generate the names file containing information about GO terms.
	# NAMESresults$GOtermNr								Numeric; GO term numbers found to be significant for input genes. 
	# NAMESresults$GOtermDescription			String; Description of GO terms = termDescription(GOtermNr). 
	# NAMESresults$GOtermId								String; GO term Id = termId(GOtermNr).
#
# Genes2GOtermsSparseMatrix		Numeric; Sparse matrix explaining the connection of genes and GO terms. 
#															Genes2GOtermsSparseMatrix[i,3]==1	iff gene in ith row is annotated to 
#															GO term in ith row.
# GO2GOAdjMatrices		List of 4:
#											Adjacency matrices for each ontology and combined sparse matrix.
	#	GO2GOAdjMatrices$GO2GOSparseAdjMatrix		Numeric; Sparse adjacency matrix describing the complete directed 
	#																					acyclic graph (DAG) of the significant GOterms up to the root, 
	#																					i.e. the edges between GOterms and their parents. Dimnames of 
	#																					columes are specifying the Ontology - BP/MF/CC.
	# GO2GOAdjMatrices$AdjMatrixGO2GOBP				Numeric; (non-sparse) Adjacency matrix of BP-DAG. 
	#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
	# GO2GOAdjMatrices$AdjMatrixGO2GOMF				Numeric; (non-sparse) Adjacency matrix of MF-DAG. 
	#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
	# GO2GOAdjMatrices$AdjMatrixGO2GOCC				Numeric; (non-sparse) Adjacency matrix of CC-DAG. 
	#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
#															
# PlotFileWithExt					String;
#													Name of the file that should be drawn with extension.
#													Extension can be one of 'png', 'eps' or 'pdf'. Default is 'png'.
#
# OPTIONAL:
# PlotDirectory						String; Default: current directory
#													Directory where PlotFileWithExt should be saved.
# MarkDetails							Boolean; Default: TRUE
#													Set TRUE if details of DAG should be marked in blue colour.
# MarkHeadlines						Boolean; Default: TRUE
#													Set TRUE if headlines should be marked in yellow colour in DAG.
# Overwrite								Boolean; Default: TRUE
#													Set TRUE if existing files with the same name should be overwritten.
#
# OUTPUT:
# One up to three DAGs each saved in seperate file.

# AUTHOR:
# CL, 20.05.2016

# USES:
# functions:
 # [1] "function"           "do.call"            "match"              "if"                 "any"               
 # [6] "is.na"              "warning"            "termId"             "length"             "print"             
# [11] "paste0"             "basename"           "file_path_sans_ext" "file_ext"           "rep"               
# [16] "as.logical"         "termDescription"    "plotGOgraph"


# Nach Ontologie aufschluesseln
AdjBP <- ORAresults$GO2GOAdjMatrices$AdjMatrixGO2GOBP
AdjMF <- ORAresults$GO2GOAdjMatrices$AdjMatrixGO2GOMF
AdjCC <- ORAresults$GO2GOAdjMatrices$AdjMatrixGO2GOCC

LrnResMatrix <- do.call(cbind, ORAresults$LRNresults)
# GOtermNr = GOtermNr, 
# OntologyNr = ontologyNr(GOtermNr), 
# NrOfGenesInUniverse = NrOfGenesInUniverse, 
# NrOfGenesInSample = NrOfGenesInSample, 
# NrOfAnnotationsInTerm = NrOfAnnotationsInTerm, 
# Up = Up, 
# ExpNrOfAnnsInTerm = ExpNrOfAnnsInTerm, 
# ObservedNrOfAnnsInTerm = ObservedNrOfAnnsInTerm, 
# RelDiff = RelDiff,
# Pvalue = Pvalues, 
# LogPvalue = log(Pvalues), 
# Certainty = Certainty, 
# InfoValue = InfoValue, 
# Remarkable = Remarkable, 
# Importance = Importance, 
# IsHeadline = IsHeadline, 
# IsDetail = IsDetail
LrnResBP <- LrnResMatrix[ORAresults$LRNresults$OntologyNr==1,,drop=FALSE]
LrnResMF <- LrnResMatrix[ORAresults$LRNresults$OntologyNr==2,,drop=FALSE]
LrnResCC <- LrnResMatrix[ORAresults$LRNresults$OntologyNr==4,,drop=FALSE]

# Finde gemeinsame GOterme in Adj und LrnRes (in Adj sind ja auch noch die ergaenzten Eltern)
BPIndex <- match(LrnResBP[,'GOtermNr'], AdjBP$GOtermNrs) # Hier duerfen keine NAs drin sein
if(any(is.na(BPIndex))){warning('drawORA: NAs gefunden. Hier ist was kaputt. Das kann nicht sein!')}
MFIndex <- match(LrnResMF[,'GOtermNr'], AdjMF$GOtermNrs) # Hier duerfen keine NAs drin sein
if(any(is.na(MFIndex))){warning('drawORA: NAs gefunden. Hier ist was kaputt. Das kann nicht sein!')}
CCIndex <- match(LrnResCC[,'GOtermNr'], AdjCC$GOtermNrs) # Hier duerfen keine NAs drin sein
if(any(is.na(CCIndex))){warning('drawORA: NAs gefunden. Hier ist was kaputt. Das kann nicht sein!')}

# Parameter fuer plotGOgraph vorbereiten:
##### BP #####
Adj <- AdjBP$AdjMatrix
GOtermIDs <- termId(AdjBP$GOtermNrs)
if(length(GOtermIDs)==0){
	print('drawORA: No GOterms in BP-DAG. Nothing to draw.')
}else{
	PlotFile <- paste0(basename(tools::file_path_sans_ext(PlotFileWithExt)), '_BP.', file_ext(PlotFileWithExt))
	# PlotDirectory <- PlotDirectory
	Significant <- rep(0,length(GOtermIDs))
	Significant[BPIndex] <- 1
	IsHeadline <- rep(FALSE, length(GOtermIDs))
	if(MarkHeadlines){
		IsHeadline[BPIndex] <- as.logical(LrnResBP[,'IsHeadline'])
	}# nur falls die markiert werden sollen den Wert angeben. Sonst alle auf Null lassen.
	# MarkDetails <- MarkDetails 
	# Overwrite <- Overwrite
	GOtermString <- termDescription(GOtermIDs)
	Remarkable <- rep(0,length(GOtermIDs))
	Remarkable[BPIndex] <- LrnResBP[,'Remarkable']
	Pvalues <- rep(1,length(GOtermIDs))
	Pvalues[BPIndex] <- LrnResBP[,'Pvalue']
	NrGenesInTerm <- rep(NaN,length(GOtermIDs))
	NrGenesInTerm[BPIndex] <- LrnResBP[,'NrOfAnnotationsInTerm']
	Expected <- rep(0,length(GOtermIDs))
	Expected[BPIndex] <- LrnResBP[,'ExpNrOfAnnsInTerm']
	Observed <- rep(0,length(GOtermIDs))
	Observed[BPIndex] <- LrnResBP[,'ObservedNrOfAnnsInTerm']
	Importance <- rep(0,length(GOtermIDs))
	Importance[BPIndex] <- LrnResBP[,'Importance']
	Up <- rep(0,length(GOtermIDs))
	Up[BPIndex] <- LrnResBP[,'Up']
	plotGOgraph(Adj, GOtermIDs, PlotFile, PlotDirectory, Significant, IsHeadline, MarkDetails, Overwrite, GOtermString, Remarkable, Pvalues,NrGenesInTerm, Expected, Observed, Importance, Up)
}# end else von if(length(GOterms)==0)

# ##### MF #####
Adj <- AdjMF$AdjMatrix
GOtermIDs <- termId(AdjMF$GOtermNrs)
if(length(GOtermIDs)==0){
	print('drawORA: No GOterms in MF-DAG. Nothing to draw.')
}else{
	PlotFile <- paste0(basename(tools::file_path_sans_ext(PlotFileWithExt)), '_MF.', file_ext(PlotFileWithExt))
	# PlotDirectory <- PlotDirectory
	Significant <- rep(0,length(GOtermIDs))
	Significant[MFIndex] <- 1
	IsHeadline <- rep(FALSE, length(GOtermIDs))
	if(MarkHeadlines){
		IsHeadline[MFIndex] <- as.logical(LrnResMF[,'IsHeadline'])
	}# nur falls die markiert werden sollen den Wert angeben. Sonst alle auf Null lassen.
	# MarkDetails <- MarkDetails 
	# Overwrite <- Overwrite
	GOtermString <- termDescription(GOtermIDs)
	Remarkable <- rep(0,length(GOtermIDs))
	Remarkable[MFIndex] <- LrnResMF[,'Remarkable']
	Pvalues <- rep(1,length(GOtermIDs))
	Pvalues[MFIndex] <- LrnResMF[,'Pvalue']
	NrGenesInTerm <- rep(NaN,length(GOtermIDs))
	NrGenesInTerm[MFIndex] <- LrnResMF[,'NrOfAnnotationsInTerm']
	Expected <- rep(0,length(GOtermIDs))
	Expected[MFIndex] <- LrnResMF[,'ExpNrOfAnnsInTerm']
	Observed <- rep(0,length(GOtermIDs))
	Observed[MFIndex] <- LrnResMF[,'ObservedNrOfAnnsInTerm']
	Importance <- rep(0,length(GOtermIDs))
	Importance[MFIndex] <- LrnResMF[,'Importance']
	Up <- rep(0,length(GOtermIDs))
	Up[MFIndex] <- LrnResMF[,'Up']
	plotGOgraph(Adj, GOtermIDs, PlotFile, PlotDirectory, Significant, IsHeadline, MarkDetails, Overwrite, GOtermString, Remarkable, Pvalues,NrGenesInTerm, Expected, Observed, Importance, Up)
}# end else von if(length(GOtermIDs)==0)

# ##### CC #####
Adj <- AdjCC$AdjMatrix
GOtermIDs <- termId(AdjCC$GOtermNrs)
if(length(GOtermIDs)==0){
	print('drawORA: No GOterms in CC-DAG. Nothing to draw.')
}else{
	PlotFile <- paste0(basename(tools::file_path_sans_ext(PlotFileWithExt)), '_CC.', file_ext(PlotFileWithExt))
	# PlotDirectory <- PlotDirectory
	Significant <- rep(0,length(GOtermIDs))
	Significant[CCIndex] <- 1
	IsHeadline <- rep(FALSE, length(GOtermIDs))
	if(MarkHeadlines){
		IsHeadline[CCIndex] <- as.logical(LrnResCC[,'IsHeadline'])
	}# nur falls die markiert werden sollen den Wert angeben. Sonst alle auf Null lassen.
	# MarkDetails <- MarkDetails 
	# Overwrite <- Overwrite
	GOtermString <- termDescription(GOtermIDs)
	Remarkable <- rep(0,length(GOtermIDs))
	Remarkable[CCIndex] <- LrnResCC[,'Remarkable']
	Pvalues <- rep(1,length(GOtermIDs))
	Pvalues[CCIndex] <- LrnResCC[,'Pvalue']
	NrGenesInTerm <- rep(NaN,length(GOtermIDs))
	NrGenesInTerm[CCIndex] <- LrnResCC[,'NrOfAnnotationsInTerm']
	Expected <- rep(0,length(GOtermIDs))
	Expected[CCIndex] <- LrnResCC[,'ExpNrOfAnnsInTerm']
	Observed <- rep(0,length(GOtermIDs))
	Observed[CCIndex] <- LrnResCC[,'ObservedNrOfAnnsInTerm']
	Importance <- rep(0,length(GOtermIDs))
	Importance[CCIndex] <- LrnResCC[,'Importance']
	Up <- rep(0,length(GOtermIDs))
	Up[CCIndex] <- LrnResCC[,'Up']
	plotGOgraph(Adj, GOtermIDs, PlotFile, PlotDirectory, Significant, IsHeadline, MarkDetails, Overwrite, GOtermString, Remarkable, Pvalues,NrGenesInTerm, Expected, Observed, Importance, Up)
}# end else von if(length(GOtermIDs)==0)

}# end function drawORA