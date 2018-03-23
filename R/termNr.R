termNr <- function(GOtermId){
# GOtermNr = termNr(GOtermId)
# die numerischen Werte der GOtermIds errechnen 
#
# INPUT
# GOtermId[1:n]             GOTermIDs,       e.g  "GO:0008150"
#
# OUTPUT
# GOtermNr[1:n]        numbers for TermIDs  e.g 8150

# ALU

#NOTA: verwendet folgende Funktionen:
#-
  
# erstmal die GO root abfangen
AllInd = which(GOtermId=="all");
if (length(AllInd)>0) { GOtermId[AllInd] <- "GO:0000000" } # all bekommt nummer 0
if(is.numeric(GOtermId)){return(GOtermId)} # Falls wir schon ne Zahl haben, gehen wir davon aus, dass das die GOtermNr ist.. 
GOtermNr = as.numeric(substr(GOtermId, 4,99));
return(GOtermNr)
} # end function termNr
