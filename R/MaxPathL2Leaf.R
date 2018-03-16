MaxPathL2Leaf <- function(OntologyName = 'BP'){
# Function to get the MaxPathL2Leaf for each term in OntologyName, 
# i.e. the length of the longest path from this node to one (arbitrary) offspring of this term.
# MaxPathL2Leaf <- MaxPathL2Leaf(OntologyName='BP')

# INPUT:
# OntologyName		Name of ontology for that the computation should be done.

# OUTPUT:
# MaxPathL2Leaf			A named vector[1:n] containing the MaxPathL2Leaf for each term in OntologyName. Names of this vector are the corresponding GO-term IDs.

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

Blaetter <- names(which(is.na(ww))) #Finde Blaetter aus der Children-Datenbank. Die, die NA sind, sind Blaetter.
MaxPathL2Leaf <- c(rep(-1,length(yy))) #Initialisiere MaxPathL2Leaf mit minus Einsen
names(MaxPathL2Leaf) <- names(yy) #MaxPathL2Leaf wird ein benannter Vektor sein, s.d. man direkt MaxPathL2Leaf und TermID zusammen hat.
count <- 0 #zaehle das MaxPathL2Leaf
MaxPathL2Leaf[Blaetter] <- 0 #Alle Blaetter-Knoten haben per Def. MaxPathL2Leaf 0
nextParentsChildren<- Blaetter
while(length(nextParentsChildren)!=0){ #solange es noch Terme gibt zu denen ich keine MaxPathL2Leaf gefunden habe machen wir weiter
	count <- count+1 #das ist dann das MaxPathL2Leaf des Knotens
	nextParentsChildren <- unique(unlist(unname(yy[nextParentsChildren])))[!is.na(unique(unlist(unname(yy[nextParentsChildren]))))] #NAs rausfiltern -> da gibts keine Eltern mehr -> also die Wurzel
	MaxPathL2Leaf[nextParentsChildren] <- count
}#end while(length(nextParentsChildren)!=0)
return(MaxPathL2Leaf = MaxPathL2Leaf)
}#end function