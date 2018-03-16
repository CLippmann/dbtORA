`ReadTXT` <-
function(FileName,InDirectory=getwd(),NrLinesToSkip=0,HasNamesInCols=FALSE,EmptySymbol='NaN'){
# V=ReadTXT(FileName,InDirectory)
# load text to a file
# INPUT
# FileName,     		name of *.txt file
# InDirectory     		InDirectory where *.txt file is
# 
# OPTIONAL
# NrLinesToSkip     so many lines are skipped bevore reading default:0
# HasNamesInCols		==TRUE means Colum 1. is a name column default FALSE
# EmptySymbol          Symbol used for NaN  default EmptySymbol='NaN'
# OUTPUT
# Data(1:n,1:d)        data i.e. all numeric colums of *.txt
# VariableNames        column names of Data
# FirstColumnNames     in case the first column consists of names these are the names, not implemented yet!

# author: MT 06/15  

CurrentDir = getwd()
setwd(InDirectory)
FileName = addext(FileName,'txt')

dataFrame = utils::read.table(FileName,skip=NrLinesToSkip,na.strings=EmptySymbol,header=HasNamesInCols)

setwd(CurrentDir)
# CL: folgende Zeile liefert bei mir nen Fehler: Hab ich deshalb auskommentiert und ohne FirstColumnNames aufgerufen.
# res=list(Data=as.matrix(dataFrame),VariableNames=colnames(as.matrix(dataFrame),FirstColumnNames=NULL))
res=list(Data=as.matrix(dataFrame),VariableNames=colnames(as.matrix(dataFrame)))

return(res)
}
