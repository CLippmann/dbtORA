WriteLRN = function(FileName, Data, Header=NULL, Key=c(), DataDefined=c(), OutDirectory=getwd(),  CommentOrDigits=NULL){
#     MATLAB! WriteLRN(FileName,Data,     Header,        Key,         Defined,        LrnDirectory, DigitsOrComment, scientific);
# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits);
# WriteLRN(FileName,Data,Header,Key,c(),OutDirectory);
# save data, column Header, and column type definition to a *.lrn file
#
# INPUT
# FileName                         name of the  file to be written
# Data(1:n,1:d)                    matrix of data , cases in rows , variables in columns, may contain NaN, NA or any other missing number symbol
#                                  for is.infinite(x) and  is.nan(x) all these x are written as NaN
#
# OPTIONAL
# OutDirectory                    the directory where to write into; if not given: current dir.
#
# Header(1:d) OR Header(1:d+1)    cell array or string matrix with column Header (d data columns) plus eventual header for Key as first entry
# Key(1:n)                        vector of row type: unique key for each line, if not given Key== [1:n]' 
# DataDefined(1:d+1)              vector of column type: keys (0=ignore,1=data,3=class,9=key)
#                                 default is leading 9 and d-times 1
#
# CommentOrDigits									either a string which is inserted as CommentOrDigits in the first line in the file
#                                 if it is a number then it is the number of significant digits (after the "."), default 6

#1.Editor: MT 3/2015 extrem wichtiger bugfix: Keyuebergabe wird nun ausgeschrieben!
# Note: Nur 13 Nachkommestellen werden ausgeschrieben, w.g write.table
# Nota: scientific not yet implemented 

checkFilename(FileName,Directory=OutDirectory,Extension='lrn',ReadOrWrite=FALSE,NameOfFunctionCalled='WriteLRN()')

if(is.list(Data)){stop('Data is a list, but it ought to be a matrix')
#WriteLRN(FileName=FileName,Data=Data$Data,Key=Data$Key,Header=Data$Header,DataDefined=Data$DataDefined,OutDirectory=OutDirectory)

}

if(is.vector(Data)){Data=as.matrix(Data)
warning('Data is a vector, but it ought to be a matrix')
}
if(!is.matrix(Data)){
stop('Data is not a matrix')
}

    
CurrentDir = getwd()
setwd(OutDirectory)

# if data is a list of length 3 -> data is the result of ReadLRN -> recursion
if (mode(as.matrix(Data))=='list' && Header(Data)==c('Data', 'Key','DataHeader', 'DataDefined', 'Key0')){
   WriteLRN(FileName, Data$Data, Data$DataNames, Data$Key, Data$DataDefined, OutDirectory)
}
else{
             

# make sure extension is there ('.lrn') 
      FileName = addext(FileName,'lrn')

#       if(is.vector(Data)){ 
#           n = length(Data)
#           m = 1
#       }else{
            n = nrow(Data)
            m = ncol(Data)
      if(missing(DataDefined)|is.null(DataDefined))
        #if(missing(Key))
             #DataDefined=c(9,rep(1,m-1))
          #else
             DataDefined=rep(1,m)
           
      if(m>=2)
            if(length(DataDefined) != m){DataDefined = c(DataDefined,rep(1,m-2))}     
#       }
      
      # ignore DataDefined==0
      
      #Data = Data[,which(DataDefined!=0)] MT: sinnlose Zeile? CL: Schmeisst man damit nicht aus Data die Spalten, die man ignorieren moechte, weil ja DataDefined==0 bedeutet, dass wir diese Daten ignorieren? Entsprechend sinnfrei finde ich die naechste Zeile. Da kuerzen wir ja nur DataDefined, s.d. es (falls wir iwo eine Null hatte) von der Laenge her nicht mehr passt und wir somit automatisch in den folgenden Zeilen mit 1 auffuellen, also keine Daten mehr ignorieren... 
      DataDefined = DataDefined[which(DataDefined!=0)]
      
      # in case DataDefined==c() or length does not fit -> only ones, i.e. defined=c(1,1,...,1)
      if(sum(DataDefined==9)>1){stop('do not indicate more than one column as key')}
      
      if(sum(DataDefined==9)==0){ 
          DataDefined = c(9,DataDefined)
      
      if(missing(Key)|is.null(Key)){#MT: Korrektur zur Keyerstellung
        Key = 1:n
        warning('Key missing, generating new Key.')
      }
      if(length(Key) != n){
			
			if(!is.null(Key)){warning('Key length is inconsistent with data length, new key is generated')}
		  }
      
          Data = cbind(as.matrix(Key), Data)
      } 
      m = ncol(Data)      
      # artificial Header if necessary

      if(length(Header)==m-1){
          aux = Header
          aux[DataDefined!=9] = Header
          aux[DataDefined==9] = 'Key'
          Header = aux
      }else if(length(Header) != m){ 
         if(!is.null(Header)){stop(paste0('Length of Header ',length(Header),' unequal length of columns ',m,'. The Length of columns is counted including the key. It is allowed to define a Header without naming the key columns with the length ',m-1))}
     
          Header= as.character(DataDefined)
          Header[DataDefined==9]='Key'
          Header[DataDefined==3]='Cls'
          Ind12 = sort(na.last=T,c(which(DataDefined==1),which(DataDefined==2)))
          Header1 = paste('C',1:length(Ind12),sep='')
          Header[Ind12] = Header1
		  }
      
      ### write in file
      # write dimensions Number of lines & columns
      header = c(paste('%\t',n),paste('%\t',length(DataDefined)))
      #MT:
      if(is.character(CommentOrDigits)){
        utils::write.table(paste0('# ',CommentOrDigits), FileName, quote=FALSE, row.names=FALSE, col.names=FALSE, na='NaN')
        utils::write.table(header, FileName,append=TRUE, quote=FALSE,sep='\t', row.names=FALSE, col.names=FALSE, na='NaN')
      }else{ 
      utils::write.table(header, FileName, quote=FALSE, row.names=FALSE, col.names=FALSE, na='NaN')
      }
      
      # write 'DataDefined'-line
      cat('% ', file=FileName, append=TRUE)
      cat(DataDefined,'\n',file=FileName, append=TRUE, sep='\t')
      
      # write 'Header'-line
      cat('% ', file=FileName, append=TRUE)
      for(i in 1:length(Header))
        Header[i]=sub(' ','',Header[i]) #Blanks ersetzen
      cat(Header,'\n', file=FileName, append=TRUE, sep='\t')
      
      if(mode(Data)!="numeric"){ #MT: Abfangen von Strings und character
        warning('Beware, non numeric values found. Trying automatic Translation into numeric data mode')
        mode(Data) <-"numeric"
      }
      # write data
      utils::write.table(Data, file=FileName, append=TRUE, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, na='NaN')

}

setwd(CurrentDir)
}
