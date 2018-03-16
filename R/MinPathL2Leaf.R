MinPathL2Leaf <- function(OntologyName = 'BP'){
# Function to get the MinPathL2Leaf for each term in OntologyName, i.e. the length 
# of the shortest path from this node to one (arbitrary) offspring of this term.
# MinPathL2Leaf <- MinPathL2Leaf(OntologyName='BP')

# INPUT:
# OntologyName		Name of ontology for that the computation should be done.

# OUTPUT:
# MinPathL2Leaf					A named vector[1:n] containing the MinPathL2Leaf for each term in OntologyName. Names of this vector are the corresponding GO-term IDs.

# USES: package GO.db

# AUTHOR:
# CL, march 2015

switch(OntologyName, #Je nach Ontologie die richtigen Datenbanken laden.
				'BP' = {ww <- as.list(GO.db::GOBPCHILDREN)
								yy <- as.list(GO.db::GOBPPARENTS)
								yy['GO:0008150']<-NA #Wir wollen 'all' nicht als oberste Wurzel haben!
								},
				'MF' = {ww <- as.list(GO.db::GOMFCHILDREN)
								yy <- as.list(GO.db::GOMFPARENTS)
								yy['GO:0003674']<-NA #Wir wollen 'all' nicht als oberste Wurzel haben!
								},
				'CC' = {ww <- as.list(GO.db::GOCCCHILDREN)
								yy <- as.list(GO.db::GOCCPARENTS)
								yy['GO:0005575']<-NA #Wir wollen 'all' nicht als oberste Wurzel haben!
								},
				stop('getLevel2 stop: OntologyName not found!')
)#end switch(OntologyName)

Blaetter <- names(which(is.na(ww))) #Finde Blaetter aus der Children-Datenbank. Die die NA sind sind Blaetter.
MinPathL2Leaf <- c(rep(-1,length(yy))) #Initialisiere MinPathL2Leaf mit minus Einsen
names(MinPathL2Leaf) <- names(yy) #MinPathL2Leaf wird ein benannter Vektor sein, s.d. man direkt MinPathL2Leaf und TermID zusammen hat.
count <- 0 #zaehle die MinPathL2Leaf
MinPathL2Leaf[Blaetter] <- 0 #Alle Blaetter-Knoten haben per Def. MinPathL2Leaf 0
nextParentsChildren<- Blaetter
while(any(MinPathL2Leaf<0)){ #solange es noch Terme gibt zu denen ich keine MinPathL2Leaf gefunden habe machen wir weiter
	count <- count+1 #das ist dann die MinPathL2Leaf des Knotens
	nextParentsChildren <- unique(unlist(unname(yy[nextParentsChildren])))[!is.na(unique(unlist(unname(yy[nextParentsChildren]))))] #NAs rausfiltern -> da gibts keine Eltern mehr -> also die Wurzel
	neue <- which(MinPathL2Leaf[nextParentsChildren] < 0) #Vorsicht: Nichts ueberschreiben!!! Die die ich schonmal gefunden habe, sollen ja die alte (erste gefundene MinPathL2Leaf) behalten.
	MinPathL2Leaf[nextParentsChildren][neue] <- count
	nextParentsChildren <- nextParentsChildren[neue]
}#end while(any(MinPathL2Leaf<0)
return(MinPathL2Leaf = MinPathL2Leaf)
}#end function