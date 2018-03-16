ReadSparseMatrix <- function(FileNameWithoutExt,InDirectory=getwd()){
# V= ReadSparseMatrix(FileNameWithoutExt,InDirectory)

# Zum Lesen von duenn besetzten Matrizen, die in *.lrn-Files gespeichert sind. (Gegebenenfalls
# mit ergaenzender *.names-Datei.)

# INPUT:
# FileNameWithoutExt						string. name of the  file to be written without extention!

# OPTIONAL:	
# InDirectory										the directory where to read from; if not given: current dir.


# OUTPUT list V with:
# ZeilenInd[1:n]								numeric vector. index of rows that contain non-zero entries
# SpaltenInd[1:n]								numeric vector. index of columns that contain non-zero entries
#																corresponding to ZeilenInd such that 
#																Matrix[ZeilenInd, SpaltenInd]==non-zero
# Inhalt[1:n]										numeric vector. Containing the values of the sparseMatrix at 
#																corresponding ZeilenInd and SpaltenInd.
# Dimension[1:2]								numeric vector. Dimensions of the complete sparseMatrix.
#																Default: c(max(ZeilenInd),max(SpaltenInd))
# Header[1:4] 									string vector. Header for columns (including Key).
# Key[1:n]											numeric vector: unique key for each line

# OPTIONAL OUTPUT: If FileNameWithoutExt.names exists:
# DimNames[1:n+x], [1:d+y]			string list. Names of the complete sparseMatrix fitting Dimension.
#																Must not be quadratic. Must not be of the same length as ZeilenInd
#																(iff there are rows containing only zeros).
#																First list element are rownames, second columnnames.
                                

# The File to read has to be a *lrn-File of the form:
# LRN named FileName in InDirectory.
# First row(s): Comment
# Following row: Number of rows and columns in FileName.lrn
# Following row: Header
# Following row:	0(Key column),
#									Number of rows in complete sparseMatrix. (Counting as well rows that do NOT contain any entry other than zero.),
# 								Number of columns in complete sparseMatrix. (Counting as well columns that do NOT contain any entry other than zero.),
#									Number of non-zero entries in complete sparseMatrix.

# If there exists a *.names-File with the same FileNameWithoutExt as the lrn-File, the *names has to be of the form:
# NAMES with rownames and columnnames in one column and a column specifiing if entry is row- or columnname.
#

# USES: ReadLRN(), possibly ReadNAMES()

# AUTHOR:
# CL 03.11.2015

Raw <- ReadLRN(paste0(FileNameWithoutExt,'.lrn'),InDirectory)
# Ueberpruefungen, ob das auch wirklich eine lrn war, in der eine sparseMatrix gespeichert war:
if(Raw$Key[1]!=0){stop('ReadSparseMatrix: Input file seems to be no sparse Matrix. Key[1] has to be 0.')}
ersteZeile <- Raw$Data[1,]
Daten <- Raw$Data[-1,]
if(dim(Daten)[2]!=3){stop('ReadSparseMatrix: Input file seems to be no sparse Matrix. Number of columns has to be 3.')}
if(max(Daten[,1]) > ersteZeile[1]){stop('ReadSparseMatrix: Input file seems to be no sparse Matrix. max(ZeilenInd) has to be smaller or equal to dimension of Matrix.')}
if(max(Daten[,2]) > ersteZeile[2]){stop('ReadSparseMatrix: Input file seems to be no sparse Matrix. max(SpaltenInd) has to be smaller or equal to dimension of Matrix.')}
if(dim(Daten)[1] != ersteZeile[3]){stop('ReadSparseMatrix: Input file seems to be no sparse Matrix. Number of non-zero entries has to be the same as there are pairs of row and column indices.')}

# So hoffen jetzt, dass es passt.
ZeilenInd <- unname(Daten[,1])
SpaltenInd <- unname(Daten[,2])
Inhalt <- unname(Daten[,3])
Dimension <- unname(ersteZeile[1:2])
Header <- c('Key', Raw$Header)
Key <- Raw$Key[-1]

# Jetzt gucken, ob wir eine *.names finden und wenn ja die Zeilen und Spaltennamen rausschreiben:
Files <- list.files(InDirectory)
if(any(Files==paste0(FileNameWithoutExt,'.names'))){
	NamesRaw <- ReadNAMES(paste0(FileNameWithoutExt,'.names'),InDirectory)
	if(is.null(NamesRaw$FurtherTexts)){# Adj-Matrix mit rownames = colnames
		DimNames <- list('rownames'=NamesRaw$Names,'colnames'=NamesRaw$Names)
	}else{
		DimNames <- list('rownames'=NamesRaw$Names[NamesRaw$FurtherTexts=='rowname'],'colnames'=NamesRaw$Names[NamesRaw$FurtherTexts=='colname'])
	}# end if(is.null(FurtherTexts))
}else{
	print('No corresponding *.names file found. No DimNames for sparseMatrix.')
	DimNames <- list('rownames'=NA, 'colnames'=NA)
}#end if

return(V=list(ZeilenInd = ZeilenInd, SpaltenInd = SpaltenInd, Inhalt = Inhalt, Dimension = Dimension, Header = Header, Key = Key, DimNames = DimNames))
}#end function