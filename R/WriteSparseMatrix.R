WriteSparseMatrix <- function(FileNameWithoutExt, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment){
# WriteSparseMatrix(FileNameWithoutExt, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, Comment){
# WriteSparseMatrix(FileNameWithoutExt, ZeilenInd, SpaltenInd, Inhalt)

# Zum Schreiben von spaerlich besetzten Matrizen. 

# Erzeugt eine *.lrn der Form:
# 	Kommentar,(ganz normal wie bei lrn)
# 	% Groesse der lrn, (ganz normal wie bei lrn)
# 	Header, (ganz normal wie bei lrn)
# 	Zeile 1 ist fuer zusaetzliche Informationen ueber Matrixdimension reserviert:
# 	In der Key-Spalte steht: 0 
# 	In der ZeilenInd-Spalte steht: Zeilenanzahl der Matrix (da es leer-Zeilen (Zeilen ohne nicht-Null-Eintraege) geben koennte)
# 	In der SpaltenInd-Spalte steht: Spaltenanzahl der Matrix (da es leer-Spalten geben koennte) 
# 	In der Inhalt-Spalte steht: Anzahl der nicht-Null-Eintraege in der gesamten Matrix.
#
# Falls DimNames gegeben wird auch noch eine *.names erzeugt, in der die Namen der gesamten Matrix
# gespeichert werden in einer Spalte. In einer zweiten Spalte wird gespeichert ob es row- oder
# columnname ist. Die Keys der lrn und der names stimmen dann selbstverstaendlich
# NICHT ueberein!!!


# INPUT
# FileNameWithoutExt						string. name of the  file to be written without extention!
# ZeilenInd[1:n]								numeric vector. index of rows that contain non-zero entries
# SpaltenInd[1:n]								numeric vector. index of columns that contain non-zero entries
#																corresponding to ZeilenInd such that 
#																Matrix[ZeilenInd, SpaltenInd]==non-zero
# Inhalt[1:n] OR [1]						numeric vector or just one number that will be repeated n-times.
#																Containing the values of the sparseMatrix at 
#																corresponding ZeilenInd and SpaltenInd.
# OPTIONAL
# Dimension[1:2]								numeric vector. If any row or column has no non-zero entries, 
#																one can specify by Dimension the original dimensions that are
#																wanted/expected.
#																Default: c(max(ZeilenInd),max(SpaltenInd))
# DimNames[1:n+x], [1:d+y]			string list. Names of the complete sparseMatrix fitting Dimension.
#																Must not be quadratic. Must not be of the same length as ZeilenInd
#																(iff there are rows containing only zeros).
#																First list element are rownames, second columnnames.
# OutDirectory									the directory where to write into; if not given: current dir.
# Header[1:4] OR [1:3]					string vector. Header for columns (including OR excluding Key).
#																Default: c('Key, ZeilenInd, SpaltenInd, Inhalt')
# Key[1:n]											numeric vector: unique key for each line, if not given Key<-1:n
# Comment												string which is inserted as comment in the first line(s) of the file
#                                

# OUTPUT:
# LRN named FileNameWithoutExt.names in OutDirectory.
# First row(s): Comment
# Following row: Number of rows and columns in FileName.lrn
# Following row: Header
# Following row:	0(Key column),
#									Number of rows in complete sparseMatrix. (Counting as well rows that do NOT contain any entry other than zero.),
# 								Number of columns in complete sparseMatrix. (Counting as well columns that do NOT contain any entry other than zero.),
#									Number of non-zero entries in complete sparseMatrix.
#
# If DimNames are specified:
# NAMES named FileNameWithoutExt.names with rownames and columnnames in one column and a column specifiing if entry is row- or columnname.

# NOTE: KEYS OF LRN AND NAMES DO NOT FIT.

# USES:
# WriteLRN(), possibly WriteNAMES()

# AUTHOR:
# CL 02.11.2015

# Ueberpruefungen:
# Laengen der uebergebenen Vektoren
if(length(ZeilenInd)!=length(SpaltenInd)){stop('WriteSparseMatrix: Numbers of entries in ZeilenInd and SpaltenInd differ but must be the same!')}#end if
if(length(Inhalt)==1){
	Inhalt <- rep(Inhalt, length(ZeilenInd))
}else{
	if(length(ZeilenInd)!=length(Inhalt)){stop('WriteSparseMatrix: Numbers of entries in ZeilenInd and Inhalt differ but must be the same or Inhalt has to be just one number!')}
}#end if
# Dimension
if(missing(Dimension)){Dimension <- c(max(ZeilenInd), max(SpaltenInd))}
if(is.null(Dimension)){Dimension <- c(max(ZeilenInd), max(SpaltenInd))}
if(length(Dimension)!=2){
	Dimension <- c(max(ZeilenInd), max(SpaltenInd))
	print('WriteSparseMatrix: Dimension has to be a numeric vector of length 2. Is now automatically set to c(max(ZeilenInd), max(SpaltenInd)).')
}#end if
# Header
if(missing(Header)){Header <- c('Key, ZeilenInd, SpaltenInd, Inhalt')}
if(is.null(Header)){Header <- c('Key, ZeilenInd, SpaltenInd, Inhalt')}
if(length(Header)==3){
	Header <- c('Key', Header)
}else{
	if(length(Header)!=4){stop('WriteSparseMatrix: Header has to be a vector of 3 or 4 entries.')}
}#end if
if(missing(Key)){Key <- 1:length(ZeilenInd)}
if(is.null(Key)){Key <- 1:length(ZeilenInd)}

ersteZeile <- c(Dimension,length(ZeilenInd))
Data <- rbind(ersteZeile, cbind(ZeilenInd, SpaltenInd, Inhalt))
KeyLong <- c(0, Key)

if(missing(OutDirectory)){OutDirectory = getwd()}
if(missing(Comment)){Comment = NULL}

# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits);
WriteLRN(paste0(FileNameWithoutExt,'.lrn'),Data,Header,KeyLong,c(),OutDirectory,Comment);

# Falls DimNames gegeben, schreiben wir auch noch eine *.names aus.
if(!missing(DimNames)){
	if(length(DimNames[[1]]) == length(DimNames[[2]]) && all(DimNames[[1]] == DimNames[[2]])){ # Falls wir eine Adj-Matrix haben, muessen wir uns nicht Spalten&Zeilennamen merken sondern nur eins.
		Names <- DimNames[[1]]
		Comments <- paste0('Names file containing names for sparse matrix saved in \n# lrn file called: ', FileNameWithoutExt, '.lrn. Row- and colnames are identical!')
		# WriteNAMES(FileName, Names,Key,FurtherTexts,OutDirectory, DescriptionHeader, Comments)
		WriteNAMES(FileNameWithoutExt, Names, 1:length(Names),OutDirectory = OutDirectory, Comments=Comments)
	}else{
		Names <- c(DimNames[[1]], DimNames[[2]])
		FurtherTexts <- c(rep('rowname',length(DimNames[[1]])),rep('colname',length(DimNames[[2]])))
		DescriptionHeader <- c('Key', 'Names', 'row/colnames')
		Comments <- paste0('Names file containing row- and columnnames for sparse matrix saved in \n# lrn file called: ', FileNameWithoutExt, '.lrn.')
		# WriteNAMES(FileName, Names,Key,FurtherTexts,OutDirectory, DescriptionHeader, Comments)
		WriteNAMES(FileNameWithoutExt, Names, 1:length(Names),FurtherTexts, OutDirectory, DescriptionHeader, Comments)
	}#end if Adjazenzmatrix
}#end if DimNames gegeben


}#end function

