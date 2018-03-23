termDescription <- function(GOtermId){
# TermStrings <- termDescription(GOtermId)
# die Strings der Terme zu den GOtermIds
#
# INPUT
# GOtermId[1:n]      GO-Term IDs,       e.g  "GO:0008150"
#
# OUTPUT
# GOTermDescription[1:n]        vector of   strings that denote the GO terms
#                           e.g. "biological process"
# EXAMPLE
# RootBPstring = termDescription("GO:0008150")   # liefert   "biological process"

#NOTA: verwendet folgende Funktionen:
# Term() aus AnnotationDbi.
#require(AnnotationDbi)
#requireNamespace(package ='AnnotationDbi', quietly = TRUE)

if(any(is.null(GOtermId)) || length(GOtermId)==0){
	if(all(is.null(GOtermId)) || length(GOtermId)==0){
		GOTermDescription <- rep(NULL, length(GOtermId))
	}else{
		NotNullInd <- which(!is.null(GOtermId))
		GOTermDescription <- rep(NULL, length(GOtermId))
		GOTermDescription[NotNullInd] <- AnnotationDbi::Term(GOtermId[NotNullInd])
	}
}else{
	GOTermDescription = AnnotationDbi::Term(GOtermId)
}
return(GOTermDescription)
}
