ReadNAMES = function(FileName,InDirectory=getwd()){
# Load *.names file
# V<- ReadNAMES(FileName,InDirectory)
# NamesKey = V$Key
# Names    = V$Names
# Rest     = V$FurtherTexts
# Comment  = V$Comments
#   
#
# INPUT
# FileName                        filename of *.names file
# OPTIONAL
# InDirectory                       InDirectory where *.names file is, default: current dir
#
# OUTPUT
# Key[1:d]                        vector, unique index from first column
# Names[1:d]                      vector, Names contained in Column 2 as char(Names), without blanks
# FurtherTexts [1:d,1:x]               vector or matrix, All in colums 3 and beyond 
# Comments [1:ccc,]               vector or matrix,  string of all lines of Comments,  without the leading "#"

#author: MT 04/2014

  if(grepl('.lrn',FileName)){
    stop('Please use ReadLRN()')
  }
  
  # if(grepl('.Data',FileName)){
  #   extensioncheck=FALSE
  #   stop('Please use ReadData()')
  # }
################################################################################################################
## Inputflow Kontrolle
################################################################################################################
  FileName  = addext(FileName,'names')
  checkFilename(FileName,Directory=InDirectory,Extension='names',ReadOrWrite=TRUE,NameOfFunctionCalled='ReadNAMES()')

  currentWD=getwd() #Aktuelles verzeichnis merken
  setwd(InDirectory) #Ins Verzeichnis wo sich Datei befindet wechseln

################################################################################################################
  #Beginn Daten einzulesen
  #Header einlesen mit HEADERSIZE
################################################################################################################
check=T
beginHeader = 1
tryCatch({
while(check){
  Header = readLines(FileName,n=beginHeader)
  if(any(grepl('%',Header))){
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

HeaderLines=readLines(FileName,n=beginHeader)
ZahlPattern = "[0-9]+"
atomicVectorInd=regexpr(ZahlPattern, HeaderLines[beginHeader])
StartRow = atomicVectorInd
EndRow = StartRow+attributes(atomicVectorInd)$match.length-1
rows=as.numeric(substr(HeaderLines[beginHeader],StartRow,EndRow))
},
error = function(c){
  warning(c)
  stop("Header or Comments are not reasonably defined, see Subversion/PUB/ZFileFormatDocuments for further instructions")
}
)

################################################################################################################
## Auslesen des Datensatzes
################################################################################################################
  #gibt fehler meldung aus wenn columns kleiner als 5 Datensätze sind.
  #Grund ist eine C-Funktion die die ersten 5 Zeilen einliest um die Daten zubestimmen
  #mehr infos: http://stackoverflow.com/questions/5990654/incomplete-final-line-warning-when-trying-to-read-a-csv-file-into-r
tryCatch({
if(!is.null(rows))  
  Data = utils::read.table(FileName,
                    sep='\t',
                    quote = "", #To disable quoting altogether
                    comment.char = "",   #turn off the interpretation of comments altogether. '#' nun auch in Daten
                    header=FALSE,         #Wenn die erste Zeile der Header ist. zB. bei csv
                    fill=TRUE, 
                    stringsAsFactors = FALSE, 
                   # na.strings=c('NaN','NA'),
                    na.strings=c(''),
                    #nrows=(rows),       #Anzahl der Zeilen zum lesen 
                    skip=beginHeader)#,    #Die ersten HEADERSIZE zeilen werden nicht nochmal gelesen
                    #colClasses=colClasses, #DatenTypen der einzlenen Spalten numeric oder character
                   # col.names=ColumnNames) #gleich die Namen der Spalten setzen
                    #colnames sollten uniqe sein, sonst wird umbenannt: mit x.1 
else
  Data = utils::read.table(FileName,
                    sep='\t',
                    quote = "", #To disable quoting altogether
                    comment.char = "",   #turn off the interpretation of comments altogether. '#' nun auch in Daten
                    header=FALSE,         #Wenn die erste Zeile der Header ist. zB. bei csv
                    fill=TRUE, 
                    stringsAsFactors = FALSE, 
                    # na.strings=c('NaN','NA'),
                    na.strings=c(''),
                    skip=beginHeader)#,    #Die ersten HEADERSIZE zeilen werden nicht nochmal gelesen
#colClasses=colClasses, #DatenTypen der einzlenen Spalten numeric oder character
# col.names=ColumnNames) #gleich die Namen der Spalten setzen
#colnames sollten uniqe sein, sonst wird umbenannt: mit x.1 

  ncols=ncol(Data)
if(!is.null(rows))
  if(rows!=nrow(Data))
    warning(paste0('ReadNAMES(): Number of rows ',nrow(Data),' does not equal number of rows in header ',rows))

################################################################################################################
## Outputflow Kontrolle
################################################################################################################

  KeyColumn=1
  Key = Data[,KeyColumn]  #Umspeichern in Key 

  Data=Data[,-(KeyColumn)]  #Key als Spalte Loeschen
  Data=data.frame(Data) #Dataframe mit einer spalte == liste -> muss wieder data.Frame sein

  #zurück wechseln ins alte Verzeichnis
  setwd(currentWD)
},
error = function(c){
  print("Key, Names or Descriuption are not reasonably defined, see Subversion/PUB/ZFileFormatDocuments for further instructions")
}
)
# Names & FurtherTexts
  NamesColumn = 2-KeyColumn
  Names=Data[,NamesColumn]
Names=as.character(Names)
Key=as.numeric(Key)
if(ncols>2){
  FurtherTexts=Data[,-NamesColumn]
  result=list(Key=Key,Names=Names,FurtherTexts=as.matrix(FurtherTexts),Comments=Comments)
}else
  result=list(Key=Key,Names=Names,Comments=Comments)

  return (result)
}
