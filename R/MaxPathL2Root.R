MaxPathL2Root <- function(OntologyName='BP'){
# to compute the MaxPathL2Root of all nodes in the three DAGs, 
# i.e. the length of the longest path from root to the term.
# MaxPathL2Root <- MaxPathL2Root(OntologyName='BP')

#INPUT:
#OntologyName			Name of the ontology we compute the MaxPathL2Root for. One of 'BP','MF' or 'CC'.

#OUTPUT:
#MaxPathL2Root		Named vector[1:n] containing the MaxPathL2Root of the GO-terms in the DAG, 
#									i.e. the maximum number of edges that connect a GO-term with the root-node. 
#									Names are corresponding GOtermIDs. Root-node has MaxPathL2Root zero.

# USES: package GO.db

#AUTHOR:
#CL, march 2015

switch(OntologyName, #Je nach Ontologie die richtige Datenbank laden und die Wurzel setzen
				'BP' = {ww <- as.list(GO.db::GOBPCHILDREN)
								Wurzel <- 'GO:0008150'
								},
				'MF' = {ww <- as.list(GO.db::GOMFCHILDREN)
								Wurzel <- 'GO:0003674'
								},
				'CC' = {ww <- as.list(GO.db::GOCCCHILDREN)
								Wurzel <- 'GO:0005575'
								},
				stop('getLevel2 stop: OntologyName not found!')
)#end switch(OntologyName)

MaxPathL2Root <- c(rep(-1,length(ww))) #MaxPathL2Root mit minus eins initialisieren
names(MaxPathL2Root) <- names(ww) #Benannten Vektor draus machen. Namen sind GOtermIDs
count <- 0 #zaehlt MaxPathL2Root
MaxPathL2Root[Wurzel] <- 0 # MaxPathL2Root der Wurzel ist per Definition Null.
nextChildsParents<- Wurzel
while(length(nextChildsParents)!=0){ #solange es noch Terme gibt zu denen ich kein MaxPathL2Root gefunden habe machen wir weiter
	count <- count+1 #das ist dann das MaxPathL2Root des Knotens
	nextChildsParents <- unique(unlist(unname(ww[nextChildsParents])))[!is.na(unique(unlist(unname(ww[nextChildsParents]))))] #Nur die, die nicht NA sind behalten. Wenn NA, dann haben wir ein Blatt gefunden.
	MaxPathL2Root[nextChildsParents] <- count
}#end while(length(nextChildsParents)!=0)
return(MaxPathL2Root = MaxPathL2Root)
}#end function