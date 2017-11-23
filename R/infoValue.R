infoValue <- function(NrOfAnnotationsInTerm, NrOfGenesInUniverse = max(NrOfAnnotationsInTerm)){
# Function calculates the partial Shannon information of gene sets in GO terms in percent,
# explaining how informative a certain term in the context of all terms is.
#
# V <- infoValue(NrOfAnnotationsInTerm, NrOfGenesInUniverse)
# InfoValue <- V$InfoValue		# information value
# InfoValueP <- V$InfoValueP	# empirical probability of occurence = NrOfAnnotationsInTerm/NrOfGenesInUniverse
#
# INPUT:
# NrOfAnnotationsInTerm			Numeric; 
#														Vector of numbers of genes annotated to corresponding GO terms.
#
# OPTIONAL:
# NrOfGenesInUniverse				Numeric; Default: max(NrOfAnnotationsInTerm)
#														Number of genes in universe. (If not restricted to reference set, NrOfGenesInUniverse
#														is the same as the number of genes (directly + indirectly) annotated to the root.)
#
# OUTPUT:
# V list of 2:
# InfoValue									Numeric;
#														A value for each term that describes how informative that term is.
# InfoValueP								Numeric;
#														The ratio of NrOfAnnotationsInTerm/NrOfGenesInUniverse for each term which is
#														the empirical probability of occurence.

# AUTHOR:
# CL, 18.03.2016 (Vorlage ALU in Matlab 2012).

# USES:
# functions:
# [1] "function" "max"      "if"       "stop"     "any"      "exp"      "log"     
# [8] "return"   "list"

# Ueberpruefe Input:
if(!is.numeric(NrOfAnnotationsInTerm)||!is.numeric(NrOfGenesInUniverse)){stop('infoValue: Input parameter have to be numeric. One ore more are not. Function stops.')}
if(length(NrOfGenesInUniverse)>1 && all(NrOfGenesInUniverse==NrOfGenesInUniverse[1])){NrOfGenesInUniverse <- NrOfGenesInUniverse[1]}
if(NrOfGenesInUniverse <= 0){stop('infoValue: NrOfGenesInUniverse has to be a positive integer greater than zero.')}
if(any(NrOfAnnotationsInTerm <= 0)){stop('infoValue: One or more values in NrOfAnnotationsInTerm are negative or zero. We are only interested in terms that have at least one (direct or indirect) annotation. Function stops.')}

InfoValueP = NrOfAnnotationsInTerm/NrOfGenesInUniverse # empirische Auftrittswkeit
InfoValue = -exp(1) * InfoValueP * log(InfoValueP) * 100 # Information im Wertebereich [0:100]

return(V = list(InfoValue = InfoValue, InfoValueP = InfoValueP))
}# end function infoValue
