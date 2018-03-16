NCBI2GeneName <- function(NCBI){
# Function to get GeneSymbol and GeneName for given NCBI numbers.

# GeneNames <- NCBI2GeneName(NCBI)

# INPUT:
# NCBI			Numeric, vector of NCBI numbers (e.g. 1).

# OUTPUT:
# list of 2
# GeneSymbol	String, vector of same length as NCBI containing abbreviation
#				of GeneName (e.g. "A1BG").
# GeneName		String, vector of same length as NCBI containingdetailed 
#				description of the gene (e.g. "alpha-1-B glycoprotein").

# AUTHOR:
# CL 27.09.2016

if(!is.numeric(NCBI)){
	stop('NCBI2GeneName: NCBI has to be a numeric vector. Function stops.')
}

Daten <- ReadNAMES('AllAnnNCBIsPlusGeneName.names', system.file('extdata',package='ORA'))
gefundeneInd <- match(NCBI, Daten$Key)
return(list(GeneSymbol = Daten$Names[gefundeneInd], GeneName = Daten$FurtherTexts[gefundeneInd]))
}# end function NCBI2GeneName