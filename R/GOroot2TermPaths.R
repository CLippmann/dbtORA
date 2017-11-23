GOroot2TermPaths <- function(TargetTerm, AdjMatrix, GOtermNr, GOroot = 8150){
# Function to get all paths from the ontology's root to one target term.

# Taxonomy <- GOroot2TermPaths(TargetTerm, AdjMatrix, GOtermNr, GOroot)

# INPUT:
# TargetTerm		Numeric;
#								GO term number of the term that the paths should be found for.
# AdjMatrix			Numeric;
#								Adjacency matrix of (only) all ancestors of TargetTerm and TargetTerm itself.
#								AdjMatrix[i,j] == 1 iff there exists an edge in 
#								GO-DAG from node i to node j (i is parent of j).
# GOtermNr			Numeric;
#								Vector of GO term numbers corresponding to rows and columns of AdjMatrix.
#								Has to contain TargetTerm.
#
# OPTIONAL:
# GOroot				Numeric; Default = 8150
#								The GO term number of the root of the ontology in which the target term 
#								and all GOtermNrs are. One of: 
#								8150 for biological process, 
#								3674 for molecular function or 
#								5575 for cellular component.
#
# OUTPUT:
# Taxonomy			Numeric;
#								A list of indices such that GOtermNr[Taxonomy[[i]]] is the i-th path from 
#								GOroot to TargetTerm.    

# AUTHOR:
# CL, 1.4. 2016, Hilfsfunktion flatten4 aus Inet   

# USES:
# functions:        
 # [1] "function"  "list"      "if"        "return"    "c"         "which"    
 # [7] "for"       "seq_along" "while"     "any"       "vapply"    "logical"  
# [13] "lapply"    "is.list"   "unlist"    "is.null"   "rownames"  "colnames" 
# [19] "names" 
# aux functions:
 #	"PfadVonNach"
 #	"flatten4"

          

# Input ueberpruefen:
if(!is.numeric(TargetTerm)||!is.numeric(GOtermNr)||!is.numeric(GOroot)||!is.numeric(AdjMatrix)){stop('GOroot2TermPaths: All input parameters have to be numeric. One or more are not! Function stops.')}
if(length(GOtermNr)!=ncol(AdjMatrix)){stop('GOroot2TermPaths: The length of vector GOtermNr has to be the same as the number of columns respectively rows of AdjMatrix. Function stops.')}
if(GOroot!=8150&&GOroot!=3674&&GOroot!=5575){stop('GOroot2TermPaths: Invalid value for GOroot. Has to be 8150 for BP, 3674 for MF or 5575 for CC. Function stops.')}

##################################
# Zwei Hilfsfunktionen definieren:

# Hilfsfunktion, um Pfade von Term "Von" zu Term "Nach" zu finden.
Pfade <- list()
PfadVonNach <- function(Von, Nach, Pfad){
if(Von == Nach){ # Abbruchbedingung! 
	return(c(Pfad, Nach)) # Wenn Von und Nach gleich sind, ist der Pfad einfach nur der eine Term.
}else{
	Pfad <- c(Pfad, Von) # An den bekannten Pfad "Von" anhaengen und zu Kindern uebergehen
	Kinder <- which(AdjMatrix[Von,] == 1) # Kinder des "Von" Terms.
	for(i in seq_along(Kinder)){ # Fuer jedes Kind des "Von" Terms, suchen wir wieder den Pfad zum "Nach" Term.
			Pfade[[i]] <- PfadVonNach(Kinder[i], Nach, Pfad) # rekursiver Aufruf der Funktion
	}# end for i
	return(Pfade) # Sehr stark verschachtelte Liste
}# end if(Von == Nach)
}# end Hilfsfunktion PfadVonNach

# Hilfsfunktion, um komplizierte Listenstruktur platt zu machen
flatten4 <- function(x){ # from http://stackoverflow.com/questions/8139677/how-to-flatten-a-list-to-a-list-without-coercion
  while(any(vapply(x, is.list, logical(1)))) { 
    x <- lapply(x, function(x) if(is.list(x)){ x }else{ list(x)})
    x <- unlist(x, recursive=FALSE) 
  }# end while
  return(x)
}# end flatten4
##################################

# Indizes des TargetTerms und der Wurzel finden:
IndexTT <- which(GOtermNr == TargetTerm)
IndexRoot <- which(GOtermNr == GOroot)

# AdjMatrix benennen, um benannte Taxonomien rauszukriegen, muessen
# dann auch die Indizes TT und Root benennen:
rownames(AdjMatrix) <- GOtermNr
colnames(AdjMatrix) <- GOtermNr
names(IndexTT) <- TargetTerm
names(IndexRoot) <- GOroot


# Pfade von Wurzel zu TargetTerm finden:
Pfad <- c()
PfadeListe <- flatten4(PfadVonNach(IndexRoot, IndexTT, Pfad))  # Beide Hilfsfunktionen nutzen

return(Taxonomy = PfadeListe)
}# end function GOroot2TermPaths
