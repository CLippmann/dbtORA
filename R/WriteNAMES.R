WriteNAMES = function(FileName, Names,  Key, FurtherTexts, OutDirectory=getwd(), DescriptionHeader, Comments){
# MATLAB: WriteNAMES(FileName,Names,Key,FurtherTexts, OutDirectory,DescriptionHeader,Comments);
   
# WriteNAMES(FileName, Names)
# WriteNAMES(FileName, Names, Key)
# WriteNAMES(FileName, Names,Key,FurtherTexts,OutDirectory, DescriptionHeader, Comments)
# Save names and eventual FurtherTexts to a *.names file.
#
# INPUT
# FileName                  name of the  file to be written
# Names[1:n]                vector with text without blanks to be put in each line
#
# OPTIONAL
# OutDirectory                 the OutDirectory where to write into; if not given: current dir.
# Key[1:n]                  vector of row type: unique key for each line, by default: [1:n]' 
# FurtherTexts[1:n,1:c]          string matrix with row FurtherTexts to be put in third column
# DescriptionHeader[1:c,]  optional header for the FurtherTexts: Array of chars, this line wil be insterted at top of file with leading #
# Comments                  vector, Char Array or matrix, these lines will be insteted at top of file with leading #
#                           Note: Also allowed: Comments='first line \n# second line'
  
#author: MT 04/2014

  
################################################################################################################
## Kontrolle des Inputflows 
################################################################################################################
keynotindata=TRUE
  colladd=1 #Anzahl die abhaengig vom Input die Anzahl an Spalten ergaenzt, abhaengig ob key als variable angeben oder in Data mit Datadefined

checkFilename(FileName,Directory=OutDirectory,Extension='names',ReadOrWrite=FALSE,NameOfFunctionCalled='WriteNAMES()')


# Pruefung welche Inputvariablen angegeben wurden
  if(!missing(Names)){ 
  Data=data.frame(Names)  
  }else{
    stop('WriteNAMES: Names ist missing')
  }
  
  if(!missing(FurtherTexts)){ 
    #MT: Anfuerungsstriche sind ein sonderzeichen, welches eine korrekte checksummen erstellen und ueberpruefung verhindert
    FurtherTexts=sub("'","",sub("'","",FurtherTexts))
    
    Data=cbind(Data,FurtherTexts)
  }

#Zeilen anzahl bestimmen
if(ncol(Data) == 0 || nrow(Data) == 0 || typeof(Data) == 'NULL'){
  stop('WriteNAMES(): No Data to write') #stop wenn Data NULL ist
}
else{
  rows = nrow(Data)
  columns = ncol(Data)
}  
    
    #richtige Datei endung
    FileName  = addext(FileName,'names')
    suppressWarnings((OutDirectory = normalizePath(OutDirectory))) 
    
    currentWD=getwd() #Aktuelles verzeichnis merken
    setwd(OutDirectory) #Ins Verzeichnis wo sich Datei befindet wechseln

################################################################################################################
## Abfang  Fehler bezueglich Key
################################################################################################################

      if(missing(Key)){ #Keyvariable nicht angegeben
        Key = 1:rows #erstellt neuen Key wenn eingebener falsch ist
        Key=as.numeric(Key)
        keynotindata=TRUE
      }else{ #Keyvariable angegeben, pruefe laenge
        if(length(Key)!=rows)stop('WriteNAMES():Length of Key isnt equal to the length of rows')
      }
 
    if(!missing(Key)){Keytmp = unique(Key)
                      if(length(Keytmp)!=length(Key)){warning(paste0('WriteNAMES(): Key with length ',length(Key),' is not unique: ',length(Keytmp)))}else{Key=Keytmp}
    } #make key uniqe
  

if(keynotindata){
    #Key hinzufuegen zu ColumnNames 
    Data = data.frame(Key,Data)
    #ColumnNames=union(KeyName,ColumnNames)
    columns=columns+1 #Wenn key Hinzugefuegt wurde, muss auch die Spalten Anzahl angepasst werden
}
    #ColumnNames zu Data
    #colnames(Data)=ColumnNames
################################################################################################################ 
    #letzter Schritt ist die Klassen mit class() anzupassen:
    #Problem war: das entweder nur string oder nur Zahl gefordert war
    #R aber keine explizite Typen definition vorraus setzt. 
################################################################################################################
    for(i in 1: columns){
      if(class(Data[,i])=='integer'){
        Data[,i]= as.numeric(Data[,i])
      }
      if(class(Data[,i])=='factor'){
        Data[,i]=as.character(Data[,i])
      }
      #else class(Data[,i])=='character' oder 'numeric' => bleibt
    }
################################################################################################################
###---------------------Beginn mit schreiben in Datei--------------------------
################################################################################################################
    #Anzahl Spalten und Zeilen schreiben
    Size=c( paste0('%\t',rows) ) 

    if(missing(Comments)){
      if(!missing(DescriptionHeader)){
        cat('#\t', file=FileName, append=TRUE)
        cat(DescriptionHeader,'\n',file=FileName, append=TRUE, sep='\t')
      }
       utils::write.table(Size, file=FileName,quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, na='NaN')

    }else{
        utils::write.table(paste0('# ',Comments), file=FileName, quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE, na='NaN')
        if(!missing(DescriptionHeader)){
          cat('#', file=FileName, append=TRUE)
          cat(DescriptionHeader,'\n',file=FileName, append=TRUE, sep='\t')
        }
        utils::write.table(Size, file=FileName,append=TRUE,quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, na='NaN')
    }


    #Daten schreiben
    utils::write.table(Data, file=FileName, append=TRUE, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, na='NaN')

##------------Daten schreiben abgeschlossen-----------  

    #zurueck wechseln ins alte Verzeichnis
    setwd(currentWD)
#  }
}
