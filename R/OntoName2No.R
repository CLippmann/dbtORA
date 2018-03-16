OntoName2No <- function(OntoName){
# Function to get the Ontology Number if input is OntologyName.
# OntoNo <- OntoName2No(OntoName)

# INPUT:
# OntoName	String vector[1:n] containing corresponding OntologyNames ('BP','MF','CC').
# OUTPUT:
# OntoNo		Numeric vector[1:n] containing Ontology Numbers (1,2 or 4).
#
# AUTHOR:
# CL 28.10.2015

# Ueberpruefe Input:
if(any(OntoName!='BP'&OntoName!='MF'&OntoName!='CC')){
	print('OntoName2No: Invalid ontology names found in input. Results in NAs in output.')
	OntoName[which(OntoName!='BP'&OntoName!='MF'&OntoName!='CC')] <- NA
}

# Wandle um:
OntoNameBP <- gsub('BP',1,OntoName)
OntoNameBPMF <- gsub('MF',2,OntoNameBP)
OntoNo <- gsub('CC',4,OntoNameBPMF)
return(as.numeric(OntoNo))
}# end function

