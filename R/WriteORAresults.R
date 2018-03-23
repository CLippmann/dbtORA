WriteORAresults <- function(FileNameWithoutExt, ORAresults, OutDirectory = getwd(), InFileWithExt = ''){
# Function to write lrn file of GO terms and the computed values, names file of GO terms and their
# descriptions and matrix of genes and GO terms where the results of the overrepresentation
# analysis (ORA) are given, the GO2GO-Matrix representing the DAG and *.names file with NCBIs
# used for computation.

# WriteORAresults(FileNameWithoutExt, ORAresults, OutDirectory)

# INPUT:
# FileNameWithoutExt		String.
#												Name of the output file without extention.
# ORAresults						List of 4. For further details see function ORA.
#												List of 18 containing results relevant for lrn:
#												LRNresults = list(GOtermNr, OntologyNr, NrOfGenesInUniverse, 
#																			NrOfGenesInSample, NrOfAnnotationsInTerm, Up, 
#																			ExpNrOfAnnsInTerm, ObservedNrOfAnnsInTerm, 
#																			RelDiff, Pvalue, LogPvalue, Certainty, InfoValue, 
#																			Remarkable, Importance, InfoContent, IsHeadline, IsDetail)
#												Names of this list will be used as Header in FileNameWithoutExt.lrn.
#												List of 3 containing results relevant for names:
#												NAMESresults = list(GOtermNr, GOtermDescription, GOtermId)
#												Names of this list will be used as Header in FileNameWithoutExt.names. 	
#												Genes2GOtermsMatrix
#												Matrix describing the annotations of genes to GO terms:
#												Dimnames will be used as Key (rownames) respectively Header (colnames).
#												List of 4 containing the adjacency matrices for each ontology and the 
#												combined sparse matrix representing the DAGs of significant GOterms:
#												GO2GOAdjMatrices = list(GO2GOSparseAdjMatrix, AdjMatrixGO2GOBP, 
#																								AdjMatrixGO2GOMF, AdjMatrixGO2GOCC)
#
# OPTIONAL:
# OutDirectory					String. Default: current directory.
#												Directory where the files should be saved.
# InFileWithExt					String. Default: ''
#												Filename of the original input file.
#
# OUTPUT:
# Files saved in OutDirectory containing entries of ORAresults.

# AUTHOR:
# CL, 8.4.2016

# USES:

# Ueberpruefe Korrektheit der Daten:
# Laengen:
if(length(ORAresults)!=4){stop('WriteORAresults: ORAresults do not have the right length. Function stops.')}
if(length(ORAresults$LRNresults)!= 19){stop('WriteORAresults: ORAresults$LRNresults are of incorrect length. Function stops.')}
if(length(ORAresults$NAMESresults)!= 3){stop('WriteORAresults: ORAresults$NAMESresults are of incorrect length. Function stops.')}
# Typen
if(!all(sapply(ORAresults$LRNresults[1:17], is.numeric))& all(sapply(ORAresults$LRNresults[18:19], is.logical)|sapply(ORAresults$LRNresults[18:19], is.numeric))){stop('WriteORAresults: One or more column(s) of ORAresults$LRNresults is not of correct type. Function stops.')}
if(!all(sapply(ORAresults$NAMESresults[2:3], mode) == 'character' &is.numeric(ORAresults$NAMESresults[[1]]))){stop('WriteORAresults: One or more column(s) of ORAresults$NAMESresults is not of correct type. Function stops.')} # sapply funktioniert nicht fuer "is.character"... -.-
SumGO2GOSparseAdjMat <- Matrix::summary(ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix)
if(!is.numeric(ORAresults$Genes2GOtermsMatrix)){stop('WriteORAresults: Entries in Genes2GOtermsMatrix must be numeric. Function stops.')}
if(!all(sapply(SumGO2GOSparseAdjMat,is.numeric))){stop('WriteORAresults: Entries in GO2GOSparseAdjMatrix must be numeric. Function stops.')}


# Jetzt ausschreiben:
# Parameter der ORA:
FilePlusParameter <- strsplit(FileNameWithoutExt, '_', TRUE)[[1]]
# Unterstriche im Dateinamen beachten:
if(FilePlusParameter[length(FilePlusParameter)] == "RefSet"){ # Wenn wir ein RefSet haben, haben wir einen Parameter mehr.
	Parameter <- paste0(tail(FilePlusParameter, n = 5), collapse = ', ')
	FileN <- paste0(head(FilePlusParameter, n = length(FilePlusParameter)-5), collapse = '_')
}else{
	Parameter <- paste0(tail(FilePlusParameter, n = 4), collapse = ', ')
	FileN <- paste0(head(FilePlusParameter, n = length(FilePlusParameter)-4), collapse = '_')
}

# Zuerst die NAMES mit allen verwendeten Genen: (NCBIs, GeneSymbols und GeneDescriptions)
ValidInputGenes <- as.numeric(dimnames(ORAresults$Genes2GOtermsMatrix)[[1]][-1])# Die GeneNCBIs sind die Zeilennamen der Genes2Terms-Matrix. In der 1.Zeile wurde eine Null ergaenzt, die muss wieder weg.
AnzValidInputGenes <- length(ValidInputGenes)
FileName <- paste0(FileN, '.names')
Names <- NCBI2GeneName(ValidInputGenes)[[2]]
Key <- ValidInputGenes
FurtherTexts <- NCBI2GeneName(ValidInputGenes)[[1]]
OutDirectory <- OutDirectory
DescriptionHeader <- c('NCBI', 'GeneSymbol', 'GeneDescription')
Comments <- paste0('Valid genes from ORA input file ', InFileWithExt, '. Duplicated or not found \n# (= nowhere annotated) genes were removed. \n# GeneDescriptions and -Symbols from "AllAnnNCBIsPlusGeneName.names" in ', system.file('extdata',package='ORA'),'.')
WriteNAMES(FileName, Names, Key, FurtherTexts, OutDirectory, DescriptionHeader, Comments)


# Die LRN mit den Termen und allen Infos:
# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)
Data <- do.call(cbind, ORAresults$LRNresults[-1]) # macht aus der Liste ne Matrix
Header <- names(ORAresults$LRNresults)
Key <- ORAresults$LRNresults[[1]]
FileName <- paste0(FileNameWithoutExt, length(Key), 'Terms.lrn')
DataDefined <- c()
OutDirectory <- OutDirectory
CommentOrDigits <- paste0('WriteORAresults: ORAresults$LRNresults. \n# ', 'Parameters: ', Parameter, '\n# ','Original input file: ', InFileWithExt)
WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)

# names zu Termen:
# WriteNAMES(FileName, Names,Key,FurtherTexts,OutDirectory, DescriptionHeader, Comments)
FileName <- paste0(FileNameWithoutExt, length(Key), 'Terms.names')
Names <- ORAresults$NAMESresults$GOtermId
NamesKey <- ORAresults$NAMESresults$GOtermNr 
FurtherTexts <- ORAresults$NAMESresults$GOtermDescription
OutDirectory <- OutDirectory
DescriptionHeader <- c('GOtermNr', 'GOtermId', 'GOtermDescription')
Comments <- paste0('WriteORAresults: ORAresults$NAMESresults. \n# ', 'Parameters: ', Parameter, '\n# ','Original input file: ', InFileWithExt)
WriteNAMES(FileName, Names,NamesKey,FurtherTexts,OutDirectory, DescriptionHeader, Comments)

if(!all(Key==NamesKey)){print('NAMES und LRN haben nicht den gleichen Key. Hier ist was kaputt!')}

# Genes2GO-Matrix
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment)
# FileName <- paste0(FileNameWithoutExt, 'GO2GenesSparseMatrix')
# ZeilenInd <- SumGenes2GOSparseMat$i
# SpaltenInd <- SumGenes2GOSparseMat$j
# Inhalt <- SumGenes2GOSparseMat$x
# Dimension <- ORAresults$Genes2GOtermsSparseMatrix@Dim
# DimNames <- ORAresults$Genes2GOtermsSparseMatrix@Dimnames
# OutDirectory <- OutDirectory
# Header <- c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt')
# Key <- seq_along(ZeilenInd)
# Comment <- 'WriteORAresults: ORAresults$Genes2GOtermsMatrix.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment)
#################################
# nicht mehr als SparseMatrix sondern als normale Matrix speichern!
# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)
Genes2GOtermsMatrix <- ORAresults$Genes2GOtermsMatrix
FileName <- paste0(FileNameWithoutExt, 'Genes2GOterms', dim(ORAresults$Genes2GOtermsMatrix)[1]-1,'x',dim(ORAresults$Genes2GOtermsMatrix)[2]) # -1 in erster Komponente, weil in der ersten Zeile ja nur die GOtermnummern stehn.
Data <- Genes2GOtermsMatrix
Header <- c('NCBIs', dimnames(Genes2GOtermsMatrix)[[2]])
Key <- as.numeric(dimnames(Genes2GOtermsMatrix)[[1]])
DataDefined <- c()
OutDirectory <- OutDirectory
CommentOrDigits <- 'WriteORAresults: GOTerms2GenesMatrix. rownames in first column: NCBI numbers, colnames in first row: GOtermNrs.'
WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)
#print(Header)

# GO2GO-Matrix
# # WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment)
# FileName <- paste0(FileNameWithoutExt, 'GO2GOSparseAdjMatrix')
# ZeilenInd <- SumGO2GOSparseAdjMat$i
# SpaltenInd <- SumGO2GOSparseAdjMat$j
# Inhalt <- SumGO2GOSparseAdjMat$x
# Dimension <- ORAresults$GO2GOAdjMatrices$GO2GOSparseAdjMatrix@Dim
# DimNames <- ORAresults$GO2GOAdjMatrices$GO2GOSparseAdjMatrix@Dimnames
# OutDirectory <- OutDirectory
# Header <- c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt')
# Key <- seq_along(ZeilenInd)
# Comment <- 'WriteORAresults: ORAresults$GO2GOSparseAdjMatrix.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment)
#################################
# nicht mehr als SparseMatrix sondern als normale Matrix speichern!
# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)
GOterms2GOtermsMatrix <- as.matrix(ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix)
FileName <- paste0(FileNameWithoutExt, 'GOterms2GOterms', ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix@Dim[1]-1,'x',ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix@Dim[2])
Data <- GOterms2GOtermsMatrix
Header <- c('GOtermNrs', dimnames(GOterms2GOtermsMatrix)[[2]])
Key <- as.numeric(dimnames(GOterms2GOtermsMatrix)[[1]])
DataDefined <- c()
OutDirectory <- OutDirectory
CommentOrDigits <- 'WriteORAresults: GOTerms2GOTerms. rownames in first column: GOtermNrs, colnames in first row: Ontology. 1==BP, 2==MF, 4==CC.'
WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits)
#print(Header)

# Passend zur GOterms-GOterms-Matrix noch ne *.names mit den GOtermen speichern:
FileName <- paste0(FileNameWithoutExt, 'GOterms', ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix@Dim[1]-1, '.names')
Key <- as.numeric(dimnames(GOterms2GOtermsMatrix)[[1]])[-1]
Names <- termId(Key)
FurtherTexts <- termDescription(Names)
OutDirectory <- OutDirectory
DescriptionHeader <- c('GOtermNr', 'GOtermID', 'GOtermDescription')
Comments <- paste0('GOterms corresponding to GOterms-GOterms-Matrix "', paste0(FileNameWithoutExt, 'GOterms2GOterms', ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix@Dim[1]-1,'x',ORAresults$GO2GOAdjMatrices$GO2GOAdjMatrix@Dim[2]), '".')
WriteNAMES(FileName, Names, Key, FurtherTexts, OutDirectory, DescriptionHeader, Comments)


}# end function WriteORAresults

