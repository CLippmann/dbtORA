drawFunctionalAreas <- function(FAfileName, FAdir, ORAfileName, AdjFileName, ORAdir = FAdir, PlotDir = FAdir, Ontology = 'BP', DrawRoot2FAs = FALSE, MarkHeadlines=TRUE, MarkDetails=TRUE, PlotFileExt='png'){
# Function draws for each Functional Area (FA) DAG with FA as root (=aspects) where ORA results 
# and FAs must be given.
# If called with missing FAfileName, filenames and directories will be requested interactively.

# drawFunctionalAreas()
# drawFunctionalAreas(Ontology = 'BP')
# drawFunctionalAreas(FAfileName, FAdir, ORAfileName, OntoADJfileName)
# drawFunctionalAreas(FAfileName, FAdir, ORAfileName, OntoADJfileName, ORAdir, PlotDir, Ontology)
# drawFunctionalAreas(FAfileName, FAdir, ORAfileName, OntoADJfileName, ORAdir, PlotDir, Ontology, DrawRoot2FAs, MarkHeadlines=TRUE, MarkDetails=TRUE, PlotFileExt='png')


# OPTIONAL INPUT:
# FAfileName				String. Filename of the *.names-File that contains the Functional Areas.
#							If not given, filenames (FAfileName, ORAfileName, OntoADJfileName) and 
#							directories (FAdir, ORAdir) will be requested interactively.
# FAdir						String. Directory where FAfileName is.
# ORAfileName				String. Filename of the *.lrn-File containing the ORA results.
# AdjFileName				String. Filename of the *.lrn-File containing the complete adjacency matrix of GOterms 
#										resultig from ORA. 
# ORAdir						String. Directory where ORAfileName and AdjFileName are.
#										Default: FAdir.
# PlotDir						String. Directory where the plots (DAGs) should be saved.
#										Default: FAdir.
#	Ontology					String. Ontology (one of: 'BP', 'CC', 'MF') for which the DAGs of FAs should be drawn.
# 									Default: 'BP'
# DrawRoot2FAs			Boolean. If DAG from Ontology root to FAs (in this Ontology) should be drawn, set TRUE.
#										Default: FALSE.
# MarkHeadlines			Boolean. If Headlines should be marked in yellow, set TRUE
#										Default: TRUE.
# MarkDetails				Boolean. If Details should be marked in blue, set TRUE.
#										Default: TRUE.
# PlotFileExt				String. Extension of the PlotFile. One of "png" or "eps".
#										Default: 'png'.


# AUTHOR: 
# CL 24.02.2016

# USES:
# package 'ggm'

# Falls keine Dateinamen angegeben wurden => interaktiv Dateipfade abfragen:
if(missing(FAfileName)){ # Wenn der erste FileName fehlt, frage alle ab.
	print('Please select *.names file containing Functional Areas.')
	flush.console() # print wird sofort ausgegeben.
	FAfile <- ask2loadFile('.names')
	FAfileName <- FAfile$FileName
	FAdir <- FAfile$InDirectory
	setwd(FAdir)
	
	print('Please select *.lrn-File containing the ORA results.')
	flush.console() # print wird sofort ausgegeben.
	ORAfile <- ask2loadFile('.lrn')
	ORAfileName <- ORAfile$FileName
	ORAdir <- ORAfile$InDirectory
	setwd(ORAdir)
	
	print('Please select *.lrn-File containing the adjacency matrix of GOterms resultig from ORA')
	print('corresponding to Ontology. Has to be in the same directory as ORA results *.lrn.')
	flush.console() # print wird sofort ausgegeben.
	AdjFileName <- ask2loadFile('.lrn')$FileName
}


#require(ggm) # fuer transClos()

# Einlesen der Daten:
# ORA:
LRN <- ReadLRN(ORAfileName , ORAdir )
# FAs:
FAreas <- ReadNAMES(FAfileName, FAdir)
# ADJ, die zu Ontology gehoert:
CompleteAdj <- ReadLRN(AdjFileName, ORAdir)
CompleteAdjZeileEins <- CompleteAdj$Data[1,] # hierdrin ist die Ontology gespeichert
CompleteAdj$Data <- CompleteAdj$Data[-1,]# hier sind nur die wirklichen Daten drin
CompleteAdj$Key <- CompleteAdj$Key[-1]
# Zeilen und Spalten aus Adj auswaehlen, die zu "Ontology" gehoeren.
Adj <- list()
switch(Ontology,
	'BP' = {Adj$Data <- CompleteAdj$Data[which(CompleteAdjZeileEins==1), which(CompleteAdjZeileEins==1)]
					Adj$Key <- CompleteAdj$Key[which(CompleteAdjZeileEins==1)]},
	'MF' = {Adj$Data <- CompleteAdj$Data[which(CompleteAdjZeileEins==2), which(CompleteAdjZeileEins==2)]
					Adj$Key <- CompleteAdj$Key[which(CompleteAdjZeileEins==2)]},
	'CC' = {Adj$Data <- CompleteAdj$Data[which(CompleteAdjZeileEins==4), which(CompleteAdjZeileEins==4)]
					Adj$Key <- CompleteAdj$Key[which(CompleteAdjZeileEins==4)]},
	stop('Ontology has to be one of "BP", "MF" or "CC".')
)
# switch(Ontology,
	# 'BP' = {Adj <- ReadLRN(BPfileADJname, ORAdir)},
	# 'MF' = {Adj <- ReadLRN(MFfileADJname, ORAdir)},
	# 'CC' = {Adj <- ReadLRN(CCfileADJname, ORAdir)},
	# stop('Ontology has to be one of "BP", "MF" or "CC".')
# )


# Bestimme die Ontologien der FAs und selektiere nur die mit "Ontology".
Ontos <- OntologyName(termId(FAreas$Key))
Keys <- FAreas$Key[Ontology==Ontos]
Names <- FAreas$Names[Ontology==Ontos]
FurtherTexts <- as.vector(FAreas$FurtherTexts)[Ontology==Ontos]

# Zeichne nur FAs. Also DAG von Wurzel bis FAs als Blaetter.
# if(DrawRoot2FAs){
# drawDAG4GOterms(GOtermId=termId(Keys), Filename='FAs', Directory = ReDi('AMLmarburg2016/20ClusterDiffEx/CL'), OutDataExtention = 'png', PrintN = TRUE, OnlyManuCur = FALSE, MakeUnique = TRUE, ColouredEdges = FALSE, EraseDOTfile = TRUE)
# }

# Zeichnen der Aspects.
FAs <- Keys
IndFAs <- match(FAs, Adj$Key)
if(any(is.na(IndFAs))){
	NAInd <- which(is.na(IndFAs))
	FAs <- FAs[-NAInd] # FAs kuerzen um die FAs, die wir nicht in der ADJ-Matrix finden koennen!
	IndFAs <- IndFAs[-NAInd] # FAs - Indizes rausschmeissen, die wir nicht finden koennen.
	warning('One or more Functional Areas were not found in AdjacencyMatrix GOtermIDs. These FAs will be ignored.')
}

Erreichbar = ggm::transClos(Adj$Data) # Erreichbar ist die transitive huelle der Adj matrix i.e. alle Pfade in Adj
# Erreichbar[i,j] = 1 <=> Von i ist j erreichbar 
# sum(Erreichbar)
# sum(Adj$Data)
# Test:
# all(colnames(Erreichbar)==termId(Adj$Key)) #TRUE :) Also nicht umsortiert. 
SpaltenInd <- list()
for(i in 1:length(FAs)){# fuer jede FA die Spalte in Erreichbar raussuchen und die Zeilenindizes der Zeilen mit 1 merken. Alle Zeilenindizes sind Nachkommen von der entsprechenden FA, zu der die Spalte gehoert.
	SpaltenInd[[i]] <- which(Erreichbar[IndFAs[i],]==1)
}
FANachkommen <- lapply(SpaltenInd, function(Si){Adj$Key[Si]}) #Term IDs der Nachkommen pro FA.
names(FANachkommen) <- FAs

# # # Mit Ergaenzung bis zur Wurzel:
# for(i in 1:length(FAs)){
# drawDAG4GOterms(GOtermId=c(termId(FANachkommen[[i]]), termId(FAs[i])), Filename=paste0('FAs',i), Directory = ReDi('AMLmarburg2016/20ClusterDiffEx/CL'), OutDataExtention = 'png', PrintN = TRUE, OnlyManuCur = FALSE, MakeUnique = TRUE, ColouredEdges = FALSE, EraseDOTfile = TRUE)
# }


# Aspects zeichen:
MarkDetails = TRUE
Overwrite = TRUE
# Filenames:
maxlength <- length(Names)
Zahlen <- sprintf(paste0("%0",max(2,nchar(maxlength)),"d"), seq_len(maxlength)) # fuelle mit fuehrenden Nullen auf
FANames <- paste0(Zahlen,strtrim(gsub("[^[:alnum:]]","",Names),20)) # nur alphanumerische Zeichen auf 20 Stellen gekuerzt
if(length(IndFAs)<length(Keys)){
	ForSequenz <- seq_along(Zahlen)[-NAInd]
}else{
	ForSequenz <- seq_along(Zahlen)
}
# GOtermIDs: (die Nachkommen und die jeweilige FA selbst pro FA als Listeneintrag)
#GOtermIDs = mapply(c, sapply(FANachkommen,termId), termId(FAs), SIMPLIFY=FALSE)
#GOtermIDs <- sapply(IndFAs, c, list(c(Adj$Key[c(IndFAs[count],SpaltenInd[[count]])])
GOtermIDs = mapply(c, termId(FAs), sapply(FANachkommen,termId), SIMPLIFY=FALSE)

# Bilde Schnitt, damit ich weiss fuer welche Terme wir Infos wie IsHeadline/Up/Significant... haben
SchnittFANachkommenUndLRN <- sapply(GOtermIDs, match, termId(LRN$Key))

	
count <- 0	
for(i in 	ForSequenz){ # Alle gefundenen FAs ausschreiben
	count <- count+1
	AdjDATA <- Adj$Data[c(IndFAs[count],SpaltenInd[[count]]),c(IndFAs[count],SpaltenInd[[count]]), drop=FALSE]# nur die zeilen und spalten von den Kindern und der FA selbst auswaehlen
	#GOtermIDs4AdjDATA <- Adj$Key[c(IndFAs[count],SpaltenInd[[count]])]
	PlotFile = paste0(FANames[i],'.', PlotFileExt)
	# jetzt aus der LRN die infos einlesen, wie ich die Knoten faerben muss
	Significant = rep(0,length(GOtermIDs[[count]])) # Initialisieren erstmal alle nicht Significant
	Significant[!is.na(SchnittFANachkommenUndLRN[[count]])] <- 1 # die, die wir im Schnitt von FANachkommen und LRN$Key gefunden haben auf Significant setzen
	IsHeadline = rep(0,length(GOtermIDs[[count]])) # erstmal alle auf "keine Headline" setzen
	IsHeadline[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],17]
	GOtermString = termDescription(GOtermIDs[[count]])
	# Remarkable = rep(-1,length(GOtermIDs[[count]])) # Remarkable auf -1 setzen. Da wo Significant == 0 wird Remarkable eh nicht ausgeschrieben!
	# Remarkable[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],14]
	Remarkable = NULL
	Pvalues = rep(-1,length(GOtermIDs[[count]])) # Pvalues auf -1 setzen.
	Pvalues[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],9]
	# NrGenesInTerm = rep(-1,length(GOtermIDs[[count]])) # NrGenesInTerm auf -1 setzen.
	# NrGenesInTerm[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],8]
	NrGenesInTerm = NULL
	Expected = rep(-1,length(GOtermIDs[[count]])) # Expected auf -1 setzen.
	Expected[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],6]
	Observed = rep(-1,length(GOtermIDs[[count]])) # Observed auf -1 setzen.
	Observed[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],7]
	Importance = NULL #rep(-1,length(GOtermIDs[[count]])) # Importance auf -1 setzen.
	# Importance[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],15]
	Up = rep(1,length(GOtermIDs[[count]])) # Alle auf 1 setzen also auf "rote Knoten" (Up wird sowieso ignoriert, wenn Significant = 0 ist!)
	Up[!is.na(SchnittFANachkommenUndLRN[[count]])] <- LRN$Data[SchnittFANachkommenUndLRN[[count]][!is.na(SchnittFANachkommenUndLRN[[count]])],5]
	# plotGOgraph(Adj,GOtermIDs,PlotFile, PlotDirectory, Significant, IsHeadline, MarkDetails, Overwrite, GOtermString, Remarkable, Pvalues, NrGenesInTerm, Expected, Observed, Importance, Up)
	plotGOgraph(AdjDATA, GOtermIDs[[count]], PlotFile, PlotDir, Significant, IsHeadline*MarkHeadlines,MarkDetails,Overwrite,GOtermString,Remarkable,Pvalues, NrGenesInTerm, Expected, Observed, Importance, Up)
}# end for plotten der gefundenen FAs

if(length(IndFAs)<length(Keys)){ # Die nicht gefundenen FAs ausschreiben. -> Bild mit nur einem Knoten mit FA-GO-TermID und Hinweis: 'Term is not a member of the GO-DAG'
	for(j in 1:length(NAInd)){
		AdjDATA <- as.matrix(0) # Adj-Matrix ist nur eine Zahl -> 0
		PlotFile = paste0(FANames[NAInd[j]],'.', PlotFileExt)
		GOtermIDs = termId(FAreas$Key[NAInd])
		Significant = 0
		IsHeadline = 0
		MarkDetails = FALSE
		Overwrite = TRUE
		GOtermString = 'Term is not a member of the GO-DAG'
		plotGOgraph(AdjDATA,GOtermIDs,PlotFile,PlotDir,Significant,IsHeadline*MarkHeadlines,MarkDetails,Overwrite,GOtermString)
	}# end for: plotte leere DAGs fuer FAs, die wir nicht gefunden haben
}# end if: falls es ungefundene FAs gibt
print('Finished.')

}#end function drawFunctionalAreas