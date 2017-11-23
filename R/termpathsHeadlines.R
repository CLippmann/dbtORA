termpathsHeadlines <- function(AdjMatrix, GOtermNr, Importance, OntologyNr = 1){
# Calculates headlines for each path from root to GO detail terms. Headlines
# represent most important nodes in these paths.
#
# PLEASE NOTE: 
# All given GO terms have to be in the specified ontology and all 
# ancestors from GO terms to ontology root (including the root-term itself) must 
# be included in the given vector of GO terms.
# 
# V <- termpathsHeadlines(AdjMatrix, GOtermNr, Importance, OntologyNr)
#
# INPUT:
# AdjMatrix			Numeric;
#								Adjacency matrix of GO terms. AdjMatrix[i,j] == 1 if there exists an edge in 
#								GO-DAG from node i to node j (i is parent of j).
# GOtermNr			Numeric;
#								Vector of all GO term numbers that are in the considered DAG. 
#								Corresponding to rows and columns of AdjMatrix.
# Importance		Numeric; Values between 0 and 100.
#								Vector of importance values (= rowwise minimum of (partial shannon) information 
#								value and certainty) for GO terms.
#
# OPTIONAL:
# OntologyNr		Numeric; Default = 1
#								To select the ontology in which the GOtermNrs are. One of: 
#								1 for biological process, 
#								2 for molecular function or 
#								4 for cellular component.
#
# OUTPUT:
# List of 3:
# Headlines					Numeric;
#										Vector of GO term numbers that are headlines, the GO term numbers with maximum
#										importance value for each path from detail to GOroot.
# AllTaxonomies			List of paths from GOroot to all details given by AdjMatrix 
#										in form of GO term numbers vectors.
#										(E.g. AllTaxonomies[[1]]==c(8150, 44699, 44763, 22402, 51231, 22) is a path from 
#										BP-root "GO:0008150" to detail "GO:0000022"via nodes 44699, 44763, 22402 and 51231
#										in BP gene ontology.)
# MaxImportanceInd	Vector of indices which indicate the position of the highest value of 
#										importance in the paths in AllTaxonomies.
#										(E.g. MaxImportanceInd==c(2,4,5) then AllTaxonomies[[1]][2] would be the GO term 
#										with the highest importance value (i.e. headline) in the first path, 
#										AllTaxonomies[[2]][4] the one in the second path and AllTaxonomies[[3]][5] the one
#										in the third path.)

# AUTHOR:
# CL, 1.4.2016

# USES:
# packages: ggm
# functions:
 # [1] "function"         "require"          "if"               "which"           
 # [5] "rowSums"          "length"           "rownames"         "colnames"        
 # [9] "lapply"           "c"                "GOroot2TermPaths" "return"          
# [13] "unlist"           "names"            "sapply"           "which.max"       
# [17] "max"              "as.numeric"       "unique"           "list"            
# [21] "unname"

# Eingabe ueberpruefen
# // TODO

#require(ggm) 
#requireNamespace(package ='ggm', quietly = TRUE)

# GOroot zur jeweiligen OntologyNr bestimmen:
if(OntologyNr == 1){GOroot <- 8150}
if(OntologyNr == 2){GOroot <- 3674}
if(OntologyNr == 4){GOroot <- 5575}

if(nrow(AdjMatrix) == 1){
	if(GOtermNr == GOroot){
	return(list(Headlines = GOtermNr, AllTaxonomies = list(GOtermNr), MaxImportanceInd = 1))
	}else{
		stop('termpathsHeadlines: AdjMatrix has dimension 1x1 but the only term is not GOroot. Please check your input. Function stops.')
	}
}

# Welche Terme sind Details:
DetailInds <- which(Matrix::rowSums(AdjMatrix) == 0)
Details <- GOtermNr[DetailInds] # alle Blaetter
AnzDetails <- length(Details)

# Welche Terme sind ueberhaupt erreichbar - ueber transitive Huelle bestimmen -> uses ggm!
Erreichbar <- ggm::transClos(AdjMatrix)

# Benenne AdjMatrix, damit man auch gleich die GOtermNrs bei den Taxos hat:
rownames(AdjMatrix) <- GOtermNr
colnames(AdjMatrix) <- GOtermNr

# Fuer jedes Detail eine eigene AdjMatrix erstellen, in der nur die von diesem Term aus 
# erreichbaren Terme bis zur Wurzel stehen. Also im Grunde eine AdjMatrix der Ancestors
# jeweils selektiert pro Detail:
PartAdjs4Details <- lapply(DetailInds, function(iterInd){AdjMatrix[c(iterInd,which(Erreichbar[,iterInd]!=0)),c(iterInd, which(Erreichbar[,iterInd]!=0)), drop = FALSE]})

# Taxonomien (Pfade) pro Detail mit Teil-AdjMatrizen berechnen:
AllTaxonomiesList <- lapply(1:AnzDetails,
											FUN=function(i,Details,GOroot,PartAdjs4Details){
												TargetTerm =  Details[i]
												Taxonomy <- GOroot2TermPaths(TargetTerm, PartAdjs4Details[[i]], as.numeric(rownames(PartAdjs4Details[[i]])), GOroot)
												return(Taxonomy)
											}, #end function in lapply
											Details, GOroot, PartAdjs4Details)#end lapply
											
# Jetzt nur noch die alle Pfade speichern - unabhaengig von Details:
AllTaxos <- lapply(unlist(AllTaxonomiesList, recursive=FALSE), names) # Liste von Vektoren mit GO term numbers

# Maximalen importance value fuer jeden Pfad finden:
names(Importance) <- GOtermNr # benennen, damit wir nachher leichter auf die richtigen values zugreifen koennen
ImpValues <- lapply(AllTaxos, function(x){Importance[x]}) # zu allen Pfaden die Importance speichern
ImpValues <- lapply(ImpValues, function(x){x[-1]}) # den ersten Eintrag jedes Pfades loeschen, da Wurzel nicht Headline sein soll
MaxImpInd <- sapply(ImpValues, which.max)
MaxImp <- sapply(ImpValues, max)
GOtermNrMaxImp <- as.numeric(names(MaxImpInd))
AllTaxos <- lapply(AllTaxos, as.numeric)

# Headlines sind diejenigen Terme mit maximaler Importance
Headlines <- unique(GOtermNrMaxImp)
return(list(Headlines = Headlines, AllTaxonomies = AllTaxos, MaxImportanceInd = unname(MaxImpInd)+1)) #+1, da bei der Berechnung die Wurzel ignoriert wurde, jetzt aber in den Pfaden wieder aufgefuehrt wird!
}# end function termpathsHeadlines
