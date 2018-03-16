ReadLRN=function(FileName=NULL,InDirectory=getwd()){
# ReadLRN reads *.lrn file   see       dbt\ZFileFormatDocuments\lrn.html
# V= ReadLRN(FileName,InDirectory)
# Key    <- V$Key     # Key[1:n]              key corresponding to the rows of Data(better be unique)
# Data   <- V$Data    # Data[1:n,1:d]         array of data: n cases in rows, d variables in columns
# Header <- V$Header  # Header[1:d,]          variable names, corresponding to the columns of Data
#
# INPUT
# FileName           name of *.lrn file
#
# OPTIONAL
# InDirectory          directory where *.lrn file is  (default ==  getwd() )
#
# OUTPUT  list with:
# Data[1:n,1:d]         array of data: n cases in rows, d variables in columns
# Key[1:n]              key corresponding to the rows of Data(better be unique)
# Header[1:d,]          variable names, corresponding to the columns of Data
# DataDefined[1:d]      the "Defined" line of *.lrn
# Comments
#
# EXAMPLE
# InDirectory  <- paste0(dbtDirectory(),'Directory');
# FileName   <- 'File.lrn';
#
# Key    <- V$Key     # Key[1:n]              key corresponding to the rows of Data(better be unique)
# Data   <- V$Data    # Data[1:n,1:d]         array of data: n cases in rows, d variables in columns
# Header <- V$Header  # Header[1:d,]          variable names, corresponding to the columns of Data
#
# author: MT 12/2015


  if(is.null(FileName)){
    res <- ask2loadFile(".lrn")
    if(is.null(res)) stop("no file selected")
    FileName = res$FileName
    InDirectory = res$InDirectory
  }  
  
FileName  = addext(FileName,'lrn');     #  checks for the right extension and adds it if necessary
# see if the diretory is there and has a canonical name
InDirectory = normalizePath(InDirectory);
checkFilename(FileName,Directory=InDirectory,Extension='lrn',ReadOrWrite=TRUE,NameOfFunctionCalled='ReadLRN()')
currentWD=getwd() #Aktuelles verzeichnis merken
setwd(InDirectory) #Ins Verzeichnis wo sich Datei befindet wechseln
Header=NULL
ColNameFlag=T
tryCatch({
check=T
beginHeader = 1
while(check){
  Header = readLines(FileName,n=beginHeader)[beginHeader]
  if(gregexpr('%',Header)[1]==1){#Prozentzeichen muss als erstes Zeichen einer Zeile kommen
    check=F
  }else{
    beginHeader = beginHeader+1 #Hier muss ich das dynamisch anpassen
  }
}

if(beginHeader>1){
  Comments=as.matrix(utils::read.table(FileName, comment.char = "%", header=FALSE,  fill=TRUE, stringsAsFactors = FALSE, na.strings=c('NA','NaN'),nrows=beginHeader-1,blank.lines.skip=T))
}else{
  Comments=NULL
}
if(!any(grepl('#',Comments))){
  Comments=NULL
}
try(
  if(!is.null(Comments)){
    Comments=sub('# ','',Comments)
    Comments=sub('#','',Comments)
    Comments=sub('\t','',Comments)
    Comments=paste(apply(Comments,1,paste,collapse=" "), collapse=" ")
  }
)

HeaderLines=readLines(FileName,n=beginHeader+4)
ZahlPattern = "[0-9]+"
atomicVectorInd=regexpr(ZahlPattern, HeaderLines[beginHeader:(beginHeader+1)])

StartRow = atomicVectorInd[1]
StartCol = atomicVectorInd[2]
EndRow = StartRow+attributes(atomicVectorInd)$match.length[1]-1
EndCol = StartCol+attributes(atomicVectorInd)$match.length[2]-1
rows=as.numeric(substr(HeaderLines[beginHeader],StartRow,EndRow))
cols=as.numeric(substr(HeaderLines[beginHeader+1],StartCol,EndCol))

DataDefinedPatternStarts = "[0-9]+"
DataDefinedInd=gregexpr(DataDefinedPatternStarts, HeaderLines[beginHeader+2])
Laengen = attributes(DataDefinedInd[[1]])$match.length
Starts =   unlist(DataDefinedInd)
Last = Starts+Laengen-1
DataDefined=as.numeric(substring(HeaderLines[beginHeader+2],Starts,Last))

Line=HeaderLines[beginHeader+3]
Line=sub('%','',Line)

HeaderPatternStarts = "([[:alnum:]]|[[:punct:]])+"
HeaderInd=gregexpr(HeaderPatternStarts,Line )
Laengen = attributes(HeaderInd[[1]])$match.length
Starts =   unlist(HeaderInd)
Last = Starts+Laengen-1
Header=substring(Line,Starts,Last)

},stop=function(f){
  warning(f)
  stop("Header or Comments are not reasonably defined, see Subversion/PUB/ZFileFormatDocuments for further instructions")
  }
)
############################
eins=length(DataDefined)
zwei=length(Header)

if(eins!=cols){
  warning(paste('Length of Datadefined',eins,'does not equal number of columns',cols))
  DataDefined=NULL
}
if(zwei!=cols){
  warning(paste('Length of Header',zwei,'does not equal number of columns',cols))
  ColNameFlag=F
}
if(eins!=zwei){
  warning(paste('Length of Datadefined',eins,'does not equal length of Header',zwei))
}
if(!is.null(DataDefined)){
  keyind=which(DataDefined==9)
}else{
    warning(paste('No Keycolumns defined! First Columns is now Keycolumn'))
    keyind=1
}  
if(length(keyind)!=1){
    warning(paste('No or More than one Keycolumns defined, number of Keykolumns: ',length(keyind)))
    warning('First Columns is now Keycolumn')
    keyind=1
}
###########################
tryCatch({
Z= scan(FileName,skip = beginHeader+3,sep='\t',what=numeric(0), quiet=T)
DataKey=matrix(Z,nrow=rows,ncol=cols,byrow=T)
Key=as.vector(DataKey[,keyind])
Data=DataKey[,-keyind]
},warning=function(e){
  warning(e)
  warning('Please check row and column numbers, Trying to bypass error...')
  Z = utils::read.table(FileName, comment.char = "%", header=FALSE,  fill=TRUE, stringsAsFactors = FALSE, na.strings=c('NA','NaN'),skip=beginHeader+3)
  Key=as.vector(Z[,keyind])
  Data=Z[,-keyind]
},stop=function(f){
  warning(f)
  stop("Header or Comments or Data are not reasonably defined, see Subversion/PUB/ZFileFormatDocuments for further instructions")
  }
)
if(!is.null(Header)){
  Header=Header[-1] #Key loeschen
  try(
  if(ColNameFlag)
    if(is.matrix(Data))
      colnames(Data)=Header
  )
}
setwd(currentWD)
res=list(Key=Key,Data=Data,Header=Header,DataDefined=DataDefined,Comments=Comments)

return(res)
}
