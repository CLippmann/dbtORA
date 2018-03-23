ontologyNr <- function(GOtermNrOrId, Verbose = FALSE){
# Returns for given GOterm (as ID or number) the corresponding gene ontology as number, where
# 1 codes for biological process, 2 for molecular function and 4 for cellular component.
# If the result is 0, something went wrong.

# OntoNr <- ontologyNr(GOtermNrOrId)
# OntoNr <- ontologyNr(GOtermNrOrId, TRUE)

# INPUT:
# GOtermNrOrId		String vector [1:n] of GOtermIds (like 'GO:0008150') OR
#									Numeric vector [1:n] of GOterm numbers (like 8150).

# OPTIONAL:
# Verbose					Boolean variable. If TRUE, function prints information in GUI window.
#									Default: FALSE

# OUTPUT:
#	OntoNumber			Numeric vector [1:n] indicating the Gene Ontology that the GOterm belongs to.
#									One of 1 for biological process, 2 for molecular function or 4 for
#									cellular component. Vector is named by GOterm IDs. Unfound GOtermNrs have
#									NaN as OntoNumber.

# USES:
# package GO.db

# AUTHOR:
#	CL, 29.06.2015


#requireNamespace(package ='GO.db', quietly = TRUE)

# Je nach Eingabe von GOterm number zu GOterm ID umwandeln (oder so lassen).
if(is.numeric(GOtermNrOrId)){
	GOtermId <- termId(GOtermNrOrId)
}else{
  GOtermId <- GOtermNrOrId					
}#end if(is.numeric(GOtermNrORId))

# Funktion Ontology() aus package GO.db aufrufen
OntoName = Ontology(GOtermId)

# Falls es welche gibt, die nicht gefunden wurden: Ausgeben (falls Verbose = TRUE ist)
NAs <- is.na(OntoName)
if(Verbose & any(is.na(OntoName))){
	warning(paste0("OntologyName: There were ", sum(NAs)," GOterm IDs that can't be found."))
	print(paste0("GOterm IDs not found: ", GOtermId[NAs]))
}# end if(Verbose & any(is.na(OntoName)))

# String in Zahl umrechnen. BP == 1, MF == 2, CC == 4.
OntoNr <- ifelse(OntoName=='BP', 1, ifelse(OntoName=='MF', 2, ifelse(OntoName=='CC', 4, NaN)))
OntoNr[NAs] <- NaN
names(OntoNr) <- GOtermId

return(OntoNumber = OntoNr)
}# end function ontologyNr