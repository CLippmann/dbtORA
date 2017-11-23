importance <- function(Certainty, InfoValue){
# Function to calculate the importance for given information values and certainty values for GO terms.

# Imp <- importance(Certainty, InfoValue)

# INPUT:
# Certainty				Numeric; Values between 0 and 100.
#									Vector of certainty values for GO terms. See certainty(Pvalues).
# InfoValue				Numeric; Values between 0 and 100.
#									Vector of (partial shannon) information values for GO terms. 
#									See infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse).
#
# OUTPUT:
# Importance			Numeric; Values between 0 and 100.
#									Vector of importance value based on InfoValue and Certainty.

# AUTHOR:
# CL, 31.3.2016

# USES:
# functions:
# [1] "function"   "if"         "is.numeric" "stop"       "any"       
# [6] "return"     "pmin" 


# Input ueberpruefen:
if(!is.numeric(InfoValue) || !is.numeric(Certainty)){stop('importance: InfoValue and Certainty have to be numeric vectors. Function stops.')}
if(any(0 > InfoValue) || any(0 > Certainty) || any(InfoValue > 100) || any(Certainty > 100)){stop('importance: InfoValue and Certainty have to be non-negative numeric vectors with values between 0 and 100. Function stops.')}
if(any(InfoValue > 100) || any(Certainty > 100)){stop('importance: Information values and/or certainty values are not in interval [0;100]. Function stops.')}

return(Importance = pmin(InfoValue, Certainty))
}# end function importance
