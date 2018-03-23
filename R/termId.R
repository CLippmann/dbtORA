termId <- function(GOtermNr){
# GOtermId <- termId(GOtermNr)
# GOtermId nach Syntax der GO fuer die nummerischen termNummern errechnen 
# GOtermId is an unique identifier for the GO term in the GO database. 
# This is always in the format GO: followed by a seven digit number 
#
# INPUT
# GOtermNr[1:n]           GO-Term NRs,          e.g. 8150 (numeric!!!)
#
# OUTPUT
# GOtermId[1:n]           GOtermIds   e.g.  "GO:0008150"

# ALU
# 1.Editors MT/CL 10/2014

#NOTA: verwendet folgende Funktionen:
#-

if(length(GOtermNr)==0){return(GOtermNr)} # Wenn man versucht von nem leeren Vektor termId aufzurufen, dann gib den leeren Vektor zurueck... 

if(is.numeric(GOtermNr[1])){
	GOtermId <- ifelse(GOtermNr<1000000, formatC(GOtermNr, digits=6, flag = 0), GOtermNr)
	GOtermId <- paste0('GO:',GOtermId)
	}else{ # not numeric, hoffen wir mal dass es schon in der endform ist
	GOtermId <- GOtermNr
}
AllInd = which(GOtermNr==0);
if (length(AllInd)>0) { GOtermId[AllInd] <- "all" } # all entspricht nummer 0

return(GOtermId);
}
