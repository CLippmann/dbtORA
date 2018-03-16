ORA <- function(NCBIs, Correction = 'BON', PvalueThreshold = 0.05, MinNrOfGenes = 2, OnlyManuCur = FALSE, RefSet = NULL, GOAall = ReadLRN('GOAall.lrn', system.file('extdata',package='ORA'))){
# Function to calculate a overrepresentation analysis, using Fisher's exact test.

# ORAResults <- ORA(NCBIs, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, RefSet)

# INPUT:
# NCBIs							Numeric;
#										Vector of NCBI numbers of genes in sample (gene set of interest).
#
# OPTIONAL:
# Correction				String; Default: 'BON'
#										Type of correction for mulitple testing of the p-values.
#										'BON' for Bonferroni, 
#										'FDR' for False Discovery Rate,
#										'RAW' if no correction should be done.
# PvalueThreshold		Numeric; Default: 0.05
#										P-value threshold. GO terms with p-values greater than PvalueThreshold will be ignored.
# MinNrOfGenes			Numeric; Default: 2
#										Minimum number of genes annotated to one term that is accepted. Only GO terms with more than
#										MinNrOfGenes  genes (of the Input genes(=NCBIs)) will be considered in calculation.
# OnlyManuCur				Boolean; Default: TRUE
#										Set TRUE if only manually curated gene annotations should be considered.
# RefSet						Numeric; Default: NULL
#										Vector of NCBI numbers of genes that form the reference set (universe). If not given i.e.
#										NULL or missing, all known genes are taken as universe.
 
# OUTPUT:
# ORAresults	list of 4:
#	LRNresults	List of 16:
#							Information needed to generate the lrn file containing all the calculated values for
#							the GO terms found to be significant for input genes.
		# LRNresults$GOtermNr								Numeric; GO term numbers found to be significant for input genes. 
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
		# LRNresults$InfoContent						Numeric; IC from GOTermInfosBP/MF/CC.lrn depending on OnlyManuCur.
		# LRNresults$InfoContentORA					Numeric; -log2(ObservedNrOfAnnsInTerm/NrOfGenesInSample)
		# LRNresults$IsHeadline							Boolean; 1 if GO term is headline, 0 if not.
		# LRNresults$IsDetail								Boolean; 1 if GO term is detail, 0 if not.
#
# NAMESresults	List of 3:
#								Information needed to generate the names file containing information about GO terms.
		# NAMESresults$GOtermNr							Numeric; GO term numbers found to be significant for input genes. 
		# NAMESresults$GOtermDescription		String; Description of GO terms = termDescription(GOtermId). 
		# NAMESresults$GOtermId							String; GO term Id = termId(GOtermNr).
#
# Genes2GOtermsMatrix									Numeric;  matrix explaining the connection of genes and GO terms. 
#														Genes2GOtermsMatrix[i,j]==1 iff gene in ith row is annotated to 
#														GO term in jth row.
# GO2GOAdjMatrices	List of 4:
#										Adjacency matrices for each ontology and combined sparse matrix.
		# GO2GOAdjMatrices$GO2GOSparseAdjMatrix		Numeric; Sparse adjacency matrix describing the complete directed 
		#																					acyclic graph (DAG) of the significant GOterms up to the root, 
		#																					i.e. the edges between GOterms and their parents. Dimnames of 
		#																					columns are specifying the ontology - BP/MF/CC.
		# GO2GOAdjMatrices$AdjMatrixGO2GOBP				Numeric; (non-sparse) Adjacency matrix of BP-DAG. 
		#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
		# GO2GOAdjMatrices$AdjMatrixGO2GOMF				Numeric; (non-sparse) Adjacency matrix of MF-DAG. 
		#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
		# GO2GOAdjMatrices$AdjMatrixGO2GOCC				Numeric; (non-sparse) Adjacency matrix of CC-DAG. 
		#																					AdjMatrixGO2GOBP[i,j]==1 iff i is parent of j.
#


# AUTHOR:
# CL 03/2016

# USES:
# packages: Matrix
# functions:
 # [1] "function"                "require"                "ReadLRN"                 "system.file"               
 # [5] "unname"                  "if"                     "which"                   "intersect"              
 # [9] "length"                  "warning"                "paste0"                  "missing"                
# [13] "is.null"                 "split"                  "lapply"                  "unique"                 
# [17] "as.character"            "unlist"                 "stop"                    "print"                  
# [21] "sapply"                  "hypergeoTest"           "exp"                     "switch"                 
# [25] "p.adjust"                "as.numeric"             "termId"                  "termNr"                 
# [29] "rep"                     "match"                  "list"                    "certainty"              
# [33] "infoValue"               "remarkableness"         "importance"              "round"                  
# [37] "termsAncestors"          "setdiff"                "union"                   "matrix"                 
# [41] "c"                       "adjMatrixTermsAncestors""rowSums"                 "names"                  
# [45] "dimnames"                "is.na"                  "termpathsHeadlines"      "return"                 
# [49] "OntologyNr"              "log"                    "termDescription" 


#require(Matrix)
#requireNamespace(package ='Matrix', quietly = TRUE)

# Finde die GO-Terme, die zu den Input-Genen gehoeren:
# Zuordnung GO-Terme zu all ihren direkten und indirekten Annotierungen aus GOdata laden:
AllNCBIs <- unname(GOAall$Data[,1])
GOtermNr <- unname(GOAall$Data[,2])
Evidence <- unname(GOAall$Data[,3])
GOsOntology <- unname(GOAall$Data[,4]) # wird nur ganz am Ende verwendet um die Ontology zu bestimmen
names(GOsOntology) <- GOtermNr

GOInfoBP <- ReadLRN('GOTermInfosBP.lrn', system.file('extdata',package='ORA'))
GOInfoMF <- ReadLRN('GOTermInfosMF.lrn', system.file('extdata',package='ORA'))
GOInfoCC <- ReadLRN('GOTermInfosCC.lrn', system.file('extdata',package='ORA'))
if(OnlyManuCur){ # IC fuer Manu und All unterschiedlich
	ICInd <- which(GOInfoBP$Header == 'InformationContentBPManu')
}else{
	ICInd <- which(GOInfoBP$Header == 'InformationContentBP')
}
GOInfoIC <- c(GOInfoBP$Data[,ICInd], GOInfoMF$Data[,ICInd], GOInfoCC$Data[,ICInd])
GOtermNrInd <- which(GOInfoBP$Header == 'GOtermNr')
GOInfoTerms <- c(GOInfoBP$Data[,GOtermNrInd], GOInfoMF$Data[,GOtermNrInd], GOInfoCC$Data[,GOtermNrInd])

if(OnlyManuCur){ #Dann brauchen wir diejenigen Terme NICHT, deren EvidenceCode "1" ist, also "IEA" lautet.
	IEATerme <- which(Evidence==1)
	AllNCBIs <- AllNCBIs[-IEATerme]
	GOtermNr <- GOtermNr[-IEATerme]
}#end if(OnlyManuCur)

ValidInputGenes <- intersect(NCBIs, AllNCBIs) # sollte optimalerweise das gleiche sein wie nur "NCBIs". Es werden nur doppelte und NCBIs, die gar nicht annotiert sind, entfernt.

# Beruecksichtige RefSet, falls gegeben, und bereite Parameter fuer hypergeometrischen Test vor:
#################################################################################################
if(missing(RefSet)|is.null(RefSet)){ # Kein RefSet gegeben:
	# Vorbereitung um Parameter fuer hypergeometischen Test zu bekommen:
	ListGenes2GOtermshlp <- split(GOtermNr, AllNCBIs) # Liste benannt mit NCBIs, beinhaltet GOTerme zu denen NCBIs annotiert sind
	ListGenes2GOterms <- lapply(ListGenes2GOtermshlp, unique)
	ListGOterms2Geneshlp <- split(AllNCBIs, GOtermNr) # Liste benannt mit GOTermen, beinhaltet annotierte NCBIs.
	ListGOterms2Genes <- lapply(ListGOterms2Geneshlp, unique)

	# Fuer Input-NCBIs diejenigen Terme finden, an denen diese annotiert sind:
	ListInputGenes2GOterms <- ListGenes2GOterms[as.character(ValidInputGenes)]
	# Zu welchen GO-Termen werden die Input-NCBIs annotiert?:
	ToInputAssociatedGOterms <- unique(unlist(unname(ListInputGenes2GOterms)))
	# Welche NCBIs aus dem Universum werden zu diesen GO-Termen (noch) annotiert?
	ListToInputAssocGOterms2Genes <- ListGOterms2Genes[as.character(ToInputAssociatedGOterms)]
	
}else{ # RefSet gegeben:
	# Vorbereitung um Parameter fuer hypergeometischen Test zu bekommen:
	IndUniverseNCBIs <- AllNCBIs %in% RefSet
	AllNCBIs <- AllNCBIs[IndUniverseNCBIs] # Aufs RefSet reduzieren
	GOtermNr <- GOtermNr[IndUniverseNCBIs] # Aufs RefSet reduzieren
	ListGenes2GOtermshlp <- split(GOtermNr, AllNCBIs) # Liste benannt mit NCBIs aus RefSet, beinhaltet GOTerme zu denen NCBIs annotiert sind
	ListGenes2GOterms <- lapply(ListGenes2GOtermshlp, unique)
	ListGOterms2Geneshlp <- split(AllNCBIs, GOtermNr) # Liste benannt mit GOTermen, beinhaltet annotierte NCBIs.
	ListGOterms2Genes <- lapply(ListGOterms2Geneshlp, unique)

	
	# Fuer Input-NCBIs diejenigen Terme finden, an denen diese annotiert sind:
	Tmp <- ValidInputGenes
	ValidInputGenes <- intersect(Tmp, RefSet)
	if(length(ValidInputGenes)==0){stop("ORA: Gene set and reference set are disjoint. No further analysis possible. Function stops.")}
	if((length(Tmp)-length(ValidInputGenes))!=0){print(paste0("ORA: ",(length(Tmp)-length(ValidInputGenes)), " gene/s is/are not in the reference set. Will be removed from input gene set."))}

	# Fuer Input-NCBIs diejenigen Terme finden, an denen diese annotiert sind:
	ListInputGenes2GOterms <- ListGenes2GOterms[as.character(ValidInputGenes)]
	# Zu welchen GO-Termen werden die Input-NCBIs annotiert?:
	ToInputAssociatedGOterms <- unique(unlist(unname(ListInputGenes2GOterms)))
	# Welche NCBIs aus dem Universum (hier RefSet) werden zu diesen GO-Termen (noch) annotiert?
	ListToInputAssocGOterms2Genes <- ListGOterms2Genes[as.character(ToInputAssociatedGOterms)]
}#end if(RefSet not given)
#################################################################################################

# Parameter fuer hypergeometrischen Test:
NrOfGenesInUniverse <- length(unique(AllNCBIs)) # Anzahl aller Gene, die zu mindestens einem GO-Term annotiert sind.
NrOfGenesInSample <- length(ValidInputGenes) # Anzahl der Gene aus dem Input minus doppelte und minus gar nicht annotierte.
NrOfAnnotationsInTerm <- sapply(ListToInputAssocGOterms2Genes, length) # Wieviele (indirekte + direkte) Annotierungen hat jeder Term (unique also ohne Evidenzen doppelt zu zaehlen!), zu dem auch mindestens ein Gen aus dem Input annotiert ist?
ObservedNrOfAnnsInTerm <- sapply(lapply(ListToInputAssocGOterms2Genes,intersect,ValidInputGenes),length) # Anzahl der Gene aus dem Input die zu den Termen annotiert sind, zu denen mindestens ein Gen aus dem Input annotiert ist.
GOtermNrs <- ToInputAssociatedGOterms

# Beruecksichtige MinNrOfGenes pro Term!
IsValidAnns <- ObservedNrOfAnnsInTerm >= MinNrOfGenes # An jedem Term im Output sollen mindestens MinNrOfGenes viele Input-Gene annotiert sein.
ObservedNrOfAnnsInTerm <- ObservedNrOfAnnsInTerm[IsValidAnns] # ObservedNrOfAnnsInTerm, die kleiner als MinNrOfGenes sind, m?ssen noch bereinigt werden
NrOfAnnotationsInTerm <- NrOfAnnotationsInTerm[IsValidAnns] # Wir brauchen auch nur die Anzahl der Annotationen fuer die Terme, fuer die gilt, dass Anzahl Observed >= Schranke ist
ExpNrOfAnnsInTerm <- NrOfAnnotationsInTerm*NrOfGenesInSample/NrOfGenesInUniverse # Erwartete Anzahl von Annotierungen pro Term ausrechnen
GOtermNrs <- GOtermNrs[IsValidAnns]
# Falls keine mehr uebrig sind nach MinNrOfGenes-Beruecksichtigung abbrechen:
if(sum(IsValidAnns)==0){stop('ORA: No GO-Terms with MinNrOfGenes annotated to them found. Function stops.')}

# Jetzt hypergeometrischen Test machen:
LogPVals <- hypergeoTest(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, NrOfGenesInUniverse, LogPvalues = TRUE)
PVals <- exp(LogPVals)


######################################################
# Korrektur fuer multiples Testen:
switch(Correction,
  'BON' = {AdjustedPvals <- p.adjust(PVals, method = 'bonferroni', n = length(PVals))},
	#Anmerkung aus Dokumentation: There seems no reason to use the unmodified Bonferroni correction 
	#because it is dominated by Holm's method, which is also valid under arbitrary assumptions. 
	'FDR' = {AdjustedPvals <- p.adjust(PVals, method = 'fdr', n = length(PVals))},
	'RAW' = {AdjustedPvals <- PVals},
	stop('ORA: Invalid input for "Correction"!')
)#end switch

######################################################
# Weitere Parameter, die zur ORA gehoeren:
Up <- as.numeric(ObservedNrOfAnnsInTerm > ExpNrOfAnnsInTerm) # if(e[i] < o[i])-> up regulated (1), else (Auch bei Gleichheit!)->  down regulated(0)

######################################################
# PvaluesThreshold beruecksichtigen:
IsRelevantPval <- AdjustedPvals <= PvalueThreshold #Boolean-Vektor welche P-Values kleiner als Pvalue-Schranke sind
GOTermId <- termId(GOtermNrs)[IsRelevantPval]  
GOtermNr <- termNr(GOTermId)
NrOfAnnotationsInTerm <- NrOfAnnotationsInTerm[IsRelevantPval]
Up <- Up[IsRelevantPval]
ExpNrOfAnnsInTerm <- ExpNrOfAnnsInTerm[IsRelevantPval]
ObservedNrOfAnnsInTerm <- ObservedNrOfAnnsInTerm[IsRelevantPval]
Pvalues <- AdjustedPvals[IsRelevantPval]
RelDiff <- (ObservedNrOfAnnsInTerm - ExpNrOfAnnsInTerm) / (0.5 * (ObservedNrOfAnnsInTerm + ExpNrOfAnnsInTerm)) * 100 #Relative Differenz in Prozent


if(length(Pvalues)==0){stop('ORA: There is no Pvalue less than the PvalueThreshold. Function stops.')}

######################################################
# output of SUMMARY in console 
print('........................................')
print('ORA: summary')
print(paste0('Number of genes in test set: ', NrOfGenesInSample))
print(paste0('Number of genes in universe/reference set: ', NrOfGenesInUniverse))
print(paste0('Number of p-values: ', length(PVals)))
print(paste0('Number of adjusted p-values: ', length(Pvalues)))


######################################################
# Matrix Genes2GOTerms erstellen:
ListInputGenes2RelevantTerms <- lapply(ListInputGenes2GOterms, intersect, GOtermNr) # Pvalue Threshold und Correction beruecksichtigen
# Schreibe als Matrix:
Laenge <- sapply(ListInputGenes2RelevantTerms, length)
gene <- rep(ValidInputGenes,Laenge)
terme <- unlist(unname(ListInputGenes2RelevantTerms))
GenIndex <- match(gene,ValidInputGenes) # ZeilenInd der sparse Matrix
TermIndex <- match(terme, GOtermNr) # SpaltenInd der sparse Matrix
#Index fuer erste Zeile - dort sollen GOtermNrs stehen:
#Zeilen <- rep(1,length(GOtermNr))
#Spalten <- seq_along(GOtermNr)
Genes2GOtermsMatrix <- matrix(0, length(ValidInputGenes)+1, length(GOtermNr))
Genes2GOtermsMatrix[cbind(GenIndex+1, TermIndex)] <- 1
Genes2GOtermsMatrix[1, ] <- GOtermNr
SpaltenNames <- gsub("[^[:alnum:]]","_",termDescription(GOTermId))
dimnames(Genes2GOtermsMatrix) <- list(c('0',ValidInputGenes), SpaltenNames)

######################################################
# Remarkableness, Certainty, Shannon Info usw. berechnen
Certainty <- certainty(Pvalues)
InfoValue <- infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse)[[1]]
Remarkable <- remarkableness(Certainty, InfoValue) # produkt Certainty & InfoValue
Importance <- importance(Certainty, InfoValue) # minimum Certainty & InfoValue
# Alle Werte in % und ohne NKS:
Certainty <- round(Certainty,2)
InfoValue <- round(InfoValue,2)
Remarkable <- round(Remarkable,2)
Importance <- round(Importance,2)

######################################################
# Adjazenzmatrix GOTerms2GOTerms des DAGs berechnen:
# Zunaechst die Vorfahren der Terme bis zur Wurzel finden und ergaenzen.
AncestorsBP <- termsAncestors(GOtermNr,1)
AncestorsMF <- termsAncestors(GOtermNr,2)
AncestorsCC <- termsAncestors(GOtermNr,4)
BPTerme <- setdiff(GOtermNr,AncestorsBP$TermsWithoutAncestors)
MFTerme <- setdiff(GOtermNr,AncestorsMF$TermsWithoutAncestors)
CCTerme <- setdiff(GOtermNr,AncestorsCC$TermsWithoutAncestors)

AllTermsBP <- union(AncestorsBP$Ancestors, BPTerme)
AllTermsMF <- union(AncestorsMF$Ancestors, MFTerme)
AllTermsCC <- union(AncestorsCC$Ancestors, CCTerme)

# Jetzt die Adj-Matrix dieser Terme suchen:
AdjMatrixGO2GOBP <- list(AdjMatrix = matrix(nrow=0, ncol=0), GOtermNrs = c())
AdjMatrixGO2GOMF <- list(AdjMatrix = matrix(nrow=0, ncol=0), GOtermNrs = c())
AdjMatrixGO2GOCC <- list(AdjMatrix = matrix(nrow=0, ncol=0), GOtermNrs = c())
BPleer <- TRUE
MFleer <- TRUE
CCleer <- TRUE
if(length(BPTerme)>0){
	AdjMatrixGO2GOBP <- adjMatrixTermsAncestors(AllTermsBP, 1)
	BPleer <- FALSE
}# end if(length(BPTerme)>0)
if(length(MFTerme)>0){
	AdjMatrixGO2GOMF <- adjMatrixTermsAncestors(AllTermsMF, 2)
	MFleer <- FALSE
}# end if(length(MFTerme)>0)
if(length(CCTerme)>0){
	AdjMatrixGO2GOCC <- adjMatrixTermsAncestors(AllTermsCC, 4)
	CCleer <- FALSE
}# end if(length(CCTerme)>0)
GO2GOSparseAdjMatrix <- Matrix::bdiag(AdjMatrixGO2GOBP$AdjMatrix, AdjMatrixGO2GOMF$AdjMatrix, AdjMatrixGO2GOCC$AdjMatrix)
ErsteZeile <- c(rep(1, dim(AdjMatrixGO2GOBP$AdjMatrix)[2]), rep(2, dim(AdjMatrixGO2GOMF$AdjMatrix)[2]), rep(4, dim(AdjMatrixGO2GOCC$AdjMatrix)[2]))
GO2GOSparseAdjMatrix <- rbind(ErsteZeile, GO2GOSparseAdjMatrix)
ZeilenDimnames <-  c('0', AdjMatrixGO2GOBP$GOtermNrs, AdjMatrixGO2GOMF$GOtermNrs, AdjMatrixGO2GOCC$GOtermNrs)
SpaltenDimnames <- strtrim(gsub("[^[:alnum:]]","",c(termDescription(termId(AdjMatrixGO2GOBP$GOtermNrs)),termDescription(termId(AdjMatrixGO2GOMF$GOtermNrs)),termDescription(termId(AdjMatrixGO2GOCC$GOtermNrs)))),20)
GO2GOSparseAdjMatrix@Dimnames <- list(ZeilenDimnames, SpaltenDimnames)

# IsDetail und IsHeadline erstellen
# Details berechnen
DetailsBP <- AdjMatrixGO2GOBP$GOtermNrs[which(Matrix::rowSums(AdjMatrixGO2GOBP$AdjMatrix)==0)]
DetailsMF <- AdjMatrixGO2GOMF$GOtermNrs[which(Matrix::rowSums(AdjMatrixGO2GOMF$AdjMatrix)==0)]
DetailsCC <- AdjMatrixGO2GOCC$GOtermNrs[which(Matrix::rowSums(AdjMatrixGO2GOCC$AdjMatrix)==0)]
AllDetails <- c(DetailsBP, DetailsMF, DetailsCC)
IsDetail <- GOtermNr %in% AllDetails



# Importance nur fuer jeweils eine Onto kuerzen und dann fuer Ancestors Nullen ergaenzen
ImpGOIDs <- termId(as.numeric(names(Importance))) # GOtermIDs der Terme, fuer die wir auch nen Pvalue haben
BPadjnames <- dimnames(AdjMatrixGO2GOBP$AdjMatrix)[[1]]
BPImpInd <- ImpGOIDs %in% BPadjnames
BPImportance <- Importance[BPImpInd]
matchindMitNA <- match(ImpGOIDs, BPadjnames)
matchindOhneNA <- matchindMitNA[!is.na(matchindMitNA)]
AllBPImportance <- rep(0,length(BPadjnames))
AllBPImportance[matchindOhneNA] <- BPImportance

ImpGOIDs <- termId(as.numeric(names(Importance))) # GOtermIDs der Terme, fuer die wir auch nen Pvalue haben
MFadjnames <- dimnames(AdjMatrixGO2GOMF$AdjMatrix)[[1]]
MFImpInd <- ImpGOIDs %in% MFadjnames
MFImportance <- Importance[MFImpInd]
matchindMitNA <- match(ImpGOIDs, MFadjnames)
matchindOhneNA <- matchindMitNA[!is.na(matchindMitNA)]
AllMFImportance <- rep(0,length(MFadjnames))
AllMFImportance[matchindOhneNA] <- MFImportance

ImpGOIDs <- termId(as.numeric(names(Importance))) # GOtermIDs der Terme, fuer die wir auch nen Pvalue haben
CCadjnames <- dimnames(AdjMatrixGO2GOCC$AdjMatrix)[[1]]
CCImpInd <- ImpGOIDs %in% CCadjnames
CCImportance <- Importance[CCImpInd]
matchindMitNA <- match(ImpGOIDs, CCadjnames)
matchindOhneNA <- matchindMitNA[!is.na(matchindMitNA)]
AllCCImportance <- rep(0,length(CCadjnames))
AllCCImportance[matchindOhneNA] <- CCImportance

# Headlines berechnen
if(!BPleer){
	TermPathsHeadlinesBP <- termpathsHeadlines(AdjMatrixGO2GOBP$AdjMatrix, AdjMatrixGO2GOBP$GOtermNrs, AllBPImportance, 1)
	BPHeadlines <- TermPathsHeadlinesBP$Headlines
}else{
	BPHeadlines <- c()
}# end if(!BPleer)
if(!MFleer){
	TermPathsHeadlinesMF <- termpathsHeadlines(AdjMatrixGO2GOMF$AdjMatrix, AdjMatrixGO2GOMF$GOtermNrs, AllMFImportance, 2)
	MFHeadlines <- TermPathsHeadlinesMF$Headlines
}else{
	MFHeadlines <- c()
}# end if(!MFleer)
if(!CCleer){
	TermPathsHeadlinesCC <- termpathsHeadlines(AdjMatrixGO2GOCC$AdjMatrix, AdjMatrixGO2GOCC$GOtermNrs, AllCCImportance, 4)
	CCHeadlines <- TermPathsHeadlinesCC$Headlines
}else{
	CCHeadlines <- c()
}# end if(!CCleer)
AllHeadlines <- c(BPHeadlines, MFHeadlines, CCHeadlines)
IsHeadline <- GOtermNr %in% AllHeadlines


Ontology <- GOsOntology[match(GOtermNr, names(GOsOntology))]

# sortieren der Eintraege: zuerst nach Ontonummer aufsteigend, dann nach Importance absteigend:
# Nur Reihenfolge merken und dann in Ausgabe ueberall hintendran haengen.
Reihenfolge <- order(Ontology,Importance*-1)

# return zusammenbauen
return(ORAresults = list(
					# Infos fuer zukuenftige LRN:
					LRNresults = list(GOtermNr = GOtermNr[Reihenfolge],  
														OntologyNr = Ontology[Reihenfolge],
														NrOfGenesInUniverse = rep(NrOfGenesInUniverse, length(Reihenfolge)), 
														NrOfGenesInSample = rep(NrOfGenesInSample, length(Reihenfolge)), 
														NrOfAnnotationsInTerm = NrOfAnnotationsInTerm[Reihenfolge], 
														Up = Up[Reihenfolge], 
														ExpNrOfAnnsInTerm = round(ExpNrOfAnnsInTerm,1)[Reihenfolge], 
														ObservedNrOfAnnsInTerm = ObservedNrOfAnnsInTerm[Reihenfolge], 
														RelDiff = round(RelDiff,1)[Reihenfolge],
														Pvalue = Pvalues[Reihenfolge], 
														LogPvalue = round(-log10(Pvalues),1)[Reihenfolge], 
														Certainty = Certainty[Reihenfolge], 
														InfoValue = InfoValue[Reihenfolge], 
														Remarkable = Remarkable[Reihenfolge], 
														Importance = Importance[Reihenfolge],
														InfoContent = GOInfoIC[match(GOtermNr,GOInfoTerms)],
														InfoContentORA = round(-log2(ObservedNrOfAnnsInTerm/NrOfGenesInSample),2),
														IsHeadline = IsHeadline[Reihenfolge], 
														IsDetail = IsDetail[Reihenfolge]), 	
					# Infos fuer zukuenftige NAMES:
					NAMESresults = list(GOtermNr = GOtermNr[Reihenfolge], 
														GOtermDescription = termDescription(termId(GOtermNr[Reihenfolge])), 
														GOtermId = termId(GOtermNr[Reihenfolge])),
					# Infos fuer Genes2GOtermsSparseMatrix:
					Genes2GOtermsMatrix = Genes2GOtermsMatrix[,Reihenfolge, drop = FALSE],
					# Infos fuer GO2GOSparseAdjMatrix:
					GO2GOAdjMatrices = list(GO2GOSparseAdjMatrix=GO2GOSparseAdjMatrix, AdjMatrixGO2GOBP=AdjMatrixGO2GOBP, AdjMatrixGO2GOMF=AdjMatrixGO2GOMF, AdjMatrixGO2GOCC=AdjMatrixGO2GOCC)
				))#end return
}# end function ORA
