hypergeoTest <- function(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, NrOfGenesInUniverse, LogPvalues = TRUE){
# Function to do a one-sided hypergeometric test.
#
# V <- hypergeoTest(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, NrOfGenesInUniverse) # Log-Pvalues
# V <- hypergeoTest(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, NrOfGenesInUniverse, LogPvalues = FALSE) # Pvalues
#
# Function to do a one-sided hypergeometric test, i.e. calculate the probability to draw more or less 
# (expectation value smaller than observed number of successes respectively expectation value greater 
# than observed number of successes) than a certain number of successes (ObservedNrOfAnnsInTerm) in a 
# fixed number of draws (NrOfGenesInSample), without replacement, from a finite population of fixed 
# size (NrOfGenesInUniverse) that contains a known number of successes (NrOfAnnotationsInTerm), wherein 
# each draw is either a success or a failure.


# INPUT:
# ObservedNrOfAnnsInTerm		Numeric;
#														Vector of observed numbers of input genes annotated to one GO term. 
# NrOfAnnotationsInTerm			Numeric;
#														Vector of numbers of all genes annotated to one GO term.
# NrOfGenesInSample					Numeric;
#														The number of input genes (genes of interest in sample) annotated to 
#														at least one GO term.
# NrOfGenesInUniverse				Numeric;
#														The number of genes in universe, i.e. all genes annotated to at 
#														least one GO term.
#
# OPTIONAL:
# LogPvalues								Boolean; Default: TRUE
#														Set TRUE if logarithmised (natural logarihm) P-values should be returned.
#
# OUTPUT:
# Pvalues										Numeric;
#														Vector of (log-)p-values of one-sided hypergeometric test. 
#														If the expected number of genes annotated to one GO term is less than 
#														ObservedNrOfAnnsInTerm, the log-p-value will be log(P(X>ObservedNrOfAnnsInTerm))
#														respectively the p-value will be P(X>ObservedNrOfAnnsInTerm),
#														where X is the hypergeometric distributed random variable. 
#														If the expected number of genes annotated to one GO term is greater than 
#														ObservedNrOfAnnsInTerm, the log-p-value will be log(P(X<ObservedNrOfAnnsInTerm)) 
#														respectively the p-value will be P(X<ObservedNrOfAnnsInTerm),
#														where X is the hypergeometric distributed random variable.
#														NOTE: In both cases the case P(X==ObservedNrOfAnnsInTerm) is excluded!

# AUTHOR:
# CL, 19.04.2016

# USES:
# functions:
 # [1] "function"         "if"               "all"              "is.numeric"      
 # [5] "c"                "stop"             "length"           "any"             
 # [9] "paste"            "which"            "LogBinomialCoeff" "rep"             
# [13] "for"              "seq_along"        "log"              "sum"             
# [17] "exp"              "return"  


# Zuerst Input ueberpruefen:
# Formal:
if(all(!is.numeric(c(ObservedNrOfAnnsInTerm, NrOfAnnotationsInTerm, NrOfGenesInSample, NrOfGenesInUniverse)))){stop('hypergeoTest: One or more input parameter are not numeric. Function stops.')}
if(length(ObservedNrOfAnnsInTerm) != length(NrOfAnnotationsInTerm)){stop('hypergeoTest: ObservedNrOfAnnsInTerm and NrOfAnnotationsInTerm need to have the same length. Function stops.')}
# Inhaltlich pruefen:
# Man darf nicht mehr Beobachtungen, als Stichproben haben: Wir koennen maximal in allen Stichproben die Beobachtung machen!
if(any(ObservedNrOfAnnsInTerm > NrOfGenesInSample)){stop('hypergeoTest: There are more observed annotations in one term than there are genes in the drawn sample. Index: ', paste(which(ObservedNrOfAnnsInTerm > NrOfGenesInSample), collapse=','),'. Function stops.')}
# Man darf nicht mehr Stichproben haben, als es Proben im Universum gibt, da wir ohne Zuruecklegen ziehen.
if(NrOfGenesInSample > NrOfGenesInUniverse){stop('hypergeoTest: There are more genes in sample than there are in universe. This is not possible. Function stops.')}
# Man darf nicht mehr Erfolge haben, als es Proben im Universum gibt.
if(any(NrOfAnnotationsInTerm > NrOfGenesInUniverse)){stop('hypergeoTest: There are more annotations in term than there are genes in universe. Index: ',paste(which(NrOfAnnotationsInTerm > NrOfGenesInUniverse), collapse=','), '. Function stops.')}
# Man kann nicht mehr Erfolge beobachten, als es insgesamt Erfolge im Universum gibt.
if(any(ObservedNrOfAnnsInTerm > NrOfAnnotationsInTerm)){stop('hypergeoTest: There are more observed annotations in one term than there are genes from universe annotated to this term. Index: ', paste(which(ObservedNrOfAnnsInTerm > NrOfAnnotationsInTerm), collapse=','), '. Function stops.')}

# Erwartungswert ausrechnen
ExpectedNrOfAnnsInTerm = (NrOfGenesInSample*NrOfAnnotationsInTerm)/NrOfGenesInUniverse

eKleinero <- ExpectedNrOfAnnsInTerm < ObservedNrOfAnnsInTerm # welche ueber-, welche unterexprimiert?
LogPvalue <- rep(-1, length(eKleinero)) # Initialisieren

# Uppertail
UpperAnn <- NrOfAnnotationsInTerm[eKleinero]
UpperObs <- ObservedNrOfAnnsInTerm[eKleinero]-1 # hier Minus eins weil wir nur P(X>x) und nicht P(X>=x) wollen
UpperGenes <- rep(NrOfGenesInSample, sum(eKleinero))
UpperUni <- rep(NrOfGenesInUniverse, sum(eKleinero))
# Lowertail
LowerAnn <- NrOfAnnotationsInTerm[!eKleinero]
LowerObs <- ObservedNrOfAnnsInTerm[!eKleinero] # hier kein Minus - es wird standardmaessig P(X<x) berechnet
LowerGenes <- rep(NrOfGenesInSample, sum(!eKleinero))
LowerUni <- rep(NrOfGenesInUniverse, sum(!eKleinero))

LogPvalue[eKleinero] <- phyper(UpperObs, UpperAnn, UpperUni-UpperAnn, UpperGenes, lower.tail = FALSE, log.p = LogPvalues)
LogPvalue[!eKleinero] <- phyper(LowerObs, LowerAnn, LowerUni-LowerAnn, LowerGenes, lower.tail = TRUE, log.p = LogPvalues)

return(Pvalues = LogPvalue)
}# end function hypergeoTest
