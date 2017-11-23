termDescription <- function(GOTermId){
# TermStrings <- termDescription(GOTermId)
# die Strings der Terme zu den GOTermIds
#
# INPUT
# GOTermId[1:n]      GO-Term IDs,       e.g  "GO:0008150"
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

if(any(is.null(GOTermId)) || length(GOTermId)==0){
	if(all(is.null(GOTermId)) || length(GOTermId)==0){
		GOTermDescription <- rep(NULL, length(GOTermId))
	}else{
		NotNullInd <- which(!is.null(GOTermId))
		GOTermDescription <- rep(NULL, length(GOTermId))
		GOTermDescription[NotNullInd] <- AnnotationDbi::Term(GOTermId[NotNullInd])
	}
}else{
	GOTermDescription = AnnotationDbi::Term(GOTermId)
}
return(GOTermDescription)
}
