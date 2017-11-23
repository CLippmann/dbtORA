termId <- function(GOTermNr){
# GOTermId <- termId(GOTermNr)
# GOTermId nach Syntax der GO fuer die nummerischen termNummern errechnen 
# GOTermId is an unique identifier for the GO term in the GO database. 
# This is always in the format GO: followed by a seven digit number 
#
# INPUT
# GOTermNr[1:n]           GO-Term NRs,          e.g. 8150 (numeric!!!)
#
# OUTPUT
# GOTermId[1:n]           GOTermIds   e.g.  "GO:0008150"

# ALU
# 1.Editors MT/CL 10/2014

#NOTA: verwendet folgende Funktionen:
#-

if(length(GOTermNr)==0){return(GOTermNr)} # Wenn man versucht von nem leeren Vektor termId aufzurufen, dann gib den leeren Vektor zurueck... 

if(is.numeric(GOTermNr[1])){
	GOTermId <- ifelse(GOTermNr<1000000, formatC(GOTermNr, digits=6, flag = 0), GOTermNr)
	GOTermId <- paste0('GO:',GOTermId)
	}else{ # not numeric, hoffen wir mal dass es schon in der endform ist
	GOTermId <- GOTermNr
}
AllInd = which(GOTermNr==0);
if (length(AllInd)>0) { GOTermId[AllInd] <- "all" } # all entspricht nummer 0

return(GOTermId);
}
