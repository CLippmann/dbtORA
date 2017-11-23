termNr <- function(GOTermId){
# GOTermNr = termNr(GOTermId)
# die numerischen Werte der GOTermIds errechnen 
#
# INPUT
# GOTermId[1:n]             GOTermIDs,       e.g  "GO:0008150"
#
# OUTPUT
# GOTermNr[1:n]        numbers for TermIDs  e.g 8150

# ALU

#NOTA: verwendet folgende Funktionen:
#-
  
# erstmal die GO root abfangen
AllInd = which(GOTermId=="all");
if (length(AllInd)>0) { GOTermId[AllInd] <- "GO:0000000" } # all bekommt nummer 0
if(is.numeric(GOTermId)){return(GOTermId)} # Falls wir schon ne Zahl haben, gehen wir davon aus, dass das die GOtermNr ist.. 
GOTermNr = as.numeric(substr(GOTermId, 4,99));
return(GOTermNr)
} # end function termNr
