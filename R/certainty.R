certainty <- function(Pvalues){
# Function to calculate the certainty (Begriffssicherheit) for a GO term, i.e. 
# the probability that there is no term with a smaller p-value in the given 
# GO subtree in percent.
#
# Certainty <- certainty(Pvalues)
#
# INPUT
# Pvalues       Numeric;
#								P-values of the GO terms.
#                         
# OUTPUT
# Certainty			Numeric;
#								Empiric probability that there is no term with a smaller p-value in the given GO subtree.

# AUTHOR:
# CL, 17.03.2016; Vorlage ALU

# USES:
# functions:
 # [1] "function"   "if"         "length"     "stop"       "all"       
 # [6] "is.numeric" "which"      "log10"      "colSums"    "sapply"    
# [11] "return"


# Input ueberpruefen:
if(length(Pvalues)==0){stop('certainty: At least one Pvalue must be given. Currently Variable "Pvalues" is empty. Function stops.')}
if(!all(is.numeric(Pvalues))){stop('certainty: Pvalues must be of type "numeric". Function stops.')}

# InformationContent ausrechnen
NichtNullInd <- which(Pvalues >= 0)
InformationContent <- Pvalues*0 # alle zuerst auf Null setzen.
# Pvalues, die nicht Null sind, in Log umrechnen fuer bessere Genauigkeit.
InformationContent[NichtNullInd] <- -log10(Pvalues[NichtNullInd])

# empirische CDF fuer InformationContent
n <- length(InformationContent)
if(n == 1){ # Falls nur ein p-Value gegeben, ist die empirische W-keit 1.
	Certainty = 1
}else{
	Certainty = Matrix::colSums(sapply(InformationContent, '>=', InformationContent))/n
}# end if(n==1)

return(Certainty*100)
}# end function certainty
