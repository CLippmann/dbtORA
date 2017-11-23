remarkableness <- function(Certainty, InfoValue){
# Function to calculate the remarkableness of a GO term.

# Rem <- remarkableness(Certainty, InfoValue)

# INPUT: 
# Certainty		Numeric;
#							(Vector of) certainty value(s) of GO term.
#							See certainty(Pvalues).
# InfoValue		Numeric;
#							(Vector of) information value(s) of GO term.
#							See infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse).
#
# OUTPUT:
# Remakable		Numeric;
#							Vector containing remarkableness values = certainty * information value of the GO term(s).

# AUTHOR:
# CL, 17.05.2016

# USES:
# functions:
# [1] "function"	"if"	"length"	"stop"	"is.numeric"	"return"

# Input ueberpruefen:
if(length(Certainty)!=length(InfoValue)){stop("remarkableness: Lengths of input parameters differ. Function stops.")}
if(!is.numeric(Certainty)&&!is.numeric(InfoValue)){stop("remarkableness: One or more input parameters are not numeric. Function stops.")}

Remarkable <- (Certainty * InfoValue)/100
return(Remarkable)
}# end function remarkableness
