MinPathL2Root <- function(OntologyName='BP'){
# to compute the MinPathL2Root of all nodes in the specified ontology, i.e. the length of the 
# shortest path form a node to the root.
# MinPathL2Root <- MinPathL2Root(OntologyName='BP')

#INPUT:
#OntologyName			Name of that ontology we compute the MinPathL2Root for. One of 'BP','MF' or 'CC'.

#OUTPUT:
#MinPathL2Root						Named vector[1:n] containing the MinPathL2Root of the GO-terms in the DAG, i.e. the minimum number of edges that connect a GO-term with the root-node. Names are corresponding GOtermIDs. Root-node has MinPathL2Root zero.

# USES: package GO.db

#AUTHOR:
#CL, march 2015

switch(OntologyName, #Je nach Ontologie die richtige Datenbank laden und Wurzel setzen
				'BP' = {ww <- as.list(GO.db::GOBPCHILDREN)
								Wurzel <- 'GO:0008150'
								},
				'MF' = {ww <- as.list(GO.db::GOMFCHILDREN)
								Wurzel <- 'GO:0003674'
								},
				'CC' = {ww <- as.list(GO.db::GOCCCHILDREN)
								Wurzel <- 'GO:0005575'
								},
				stop('getDepth2 stop: OntologyName not found!')
)#end switch(OntologyName)

MinPathL2Root <- c(rep(-1,length(ww))) #Initialisiere mit minus Eins
names(MinPathL2Root) <- names(ww) #Benannter Vektor.
count <- 0 #zaehlt MinPathL2Root
MinPathL2Root[Wurzel] <- 0 # MinPathL2Root der Wurzel per Def. = 0.
nextChildsParents<- Wurzel
while(any(MinPathL2Root<0)){ #solange es noch Terme gibt zu denen ich keine MinPathL2Root gefunden habe machen wir weiter
	count <- count+1 #das ist dann die MinPathL2Root des Knotens
	nextChildsParents <- unique(unlist(unname(ww[nextChildsParents])))[!is.na(unique(unlist(unname(ww[nextChildsParents]))))] # Nur die, die nicht NA sind nehmen. Alle Blaetter haben NAs.
	neue <- which(MinPathL2Root[nextChildsParents] < 0) #Vorsicht: nicht ueberschreiben.
	MinPathL2Root[nextChildsParents][neue] <- count
	nextChildsParents <- nextChildsParents[neue]
}#end while(any(MinPathL2Root<0)
return(MinPathL2Root = MinPathL2Root)
}#end function