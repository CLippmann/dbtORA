plotGOgraph <- function(Adj,GOtermIDs,PlotFile,PlotDirectory=getwd(),
													Significant=rep(1,length(GOtermIDs)),IsHeadline=rep(0,length(GOtermIDs)),
													MarkDetails=TRUE, Overwrite=TRUE, GOtermString=NULL,Remarkable=NULL,Pvalues=NULL,
													NrGenesInTerm=NULL,Expected = NULL, Observed = NULL, Importance = NULL, Up=NULL){
# plotGOgraph(Adj,GOtermIDs,PlotFile,PlotDirectory,Significant,IsHeadline,MarkDetails,Overwrite,GOtermString,Remarkable,Pvalues,NrGenesInTerm,Up)
# plotGOgraph(Adj,GOtermIDs,PlotFile)

# zeichnen des GO Graphen und Ausgabe der Zeichnung in eine Datei (Dateiformat durch Endung in 
# PlotFile bestimmt).
#
# INPUT:
# Adj[1:n,1:n]						Adjazenzmatrix  passend zu den folgenden Knoten Infos mit Adj[i,j]==1 
#													iff i is parent of j. Muss ganzen DAG enthalten, damit ein sinnvolles Ergebnis
#													entsteht.
# GOtermIDs[1:n]					GO-Term IDs,     e.g  "GO:0008150", inklusive Eltern bis zur Wurzel.
# PlotFile								Name der Datei in die ausgegeben wird mit Dateiendung. 
#													Falls keine Dateiendung gegeben ist, wird ein png ausgegeben.
#
# OPTIONAL:
# PlotDirectory						Ausgabedirectory. default = getwd()    
# Significant[1:n]				Knoten mit Significant==1 werden rot, die anderen weiss gezeichntet 
#													(default == ones())
# IsHeadline[1:n]					Knoten mit IsHeadline==1 werden gelb markiert (default == 0)
# MarkDetails							TRUE (default) means Details = Blaetter werden blau gezeichnet
# Overwrite								TRUE, falls Files in PlotDirectory mit dem gleichen Dateinamen ueberschrieben 
#													werden sollen. Falls FALSE, wird fortlaufende Nummer angefuegt. Default: TRUE
# GOtermString[1:n]				GO-Term strings, e.g   "biological_process". Falls nicht gegeben (NULL), 
#													werden GOtermStrings automatisch bestimmt.
# Remarkable[1:n]					die Zahl, welche im Knoten gezeigt wird. Falls nicht gegeben(NULL), wird 
#													keine Zahl im Knoten angezeigt.
# Pvalues[1:n]						falls angegeben, wird auch dieser Wert in die Knoten gezeichnet. Default: NULL
# NrGenesInTerm[1:n]			falls angegeben, wird auch dieser Wert in die Knoten gezeichnet. Default: NULL
# Observed[1:n]						Anzahl beobachteter Annotierungen am Term, falls gegeben wird der Wert auch in 
#													den Knoten geschrieben. Default: NULL.
# Expected[1:n]						Anzahl statistisch erwarteter Annotierungen am Term. Falls gegeben, wird der Wert
#													auch in den Knoten geschrieben. Default: NULL.
#	Importance[1:n]					falls angegeben, wird auch dieser Wert in die Knoten gezeichnet. Default: NULL.
# Up[1:n]									Are there more than the expected nr of genes in the GOterm, dann 1 sonst 0.
# 												Knoten mit which(Up==1) werden rot, die andereen gruen gezeichnet. Default: NULL

# USES:
# package GO.db
# functions termDescription, checkFilename, fileparts

# AUTHOR:
# CL
# Edit: Man muss keine Zahlen (Remarkable, Pvalues etc.) mehr angeben sondern kann auch einfach 
#				nur Farben Text und GOtermID zeichnen.

# Parameter ueberpruefen:

# Da die Adj frueher immer andersrum war, einfach hier einmal umdrehen um nicht alles umschreiben zu muessen:
Adj <- t(Adj)

# Adj quadratisch? 
if(nrow(Adj) != ncol(Adj)){stop('plotGOgraph: Adjacency matrix has to be quadratic! Function stops.')}
# GOtermIDs genauso lang wie Zeilen- bzw. Spaltenanzahl in Adj?
if(nrow(Adj) != length(GOtermIDs)){stop('plotGOgraph: Length of GOtermIDs differs from row and column dimension of adjacency matrix. Function stops.')}
# GOtermIDs immer 10 character lang und tatsaechlich IDs und keine Nummern?
if(any(nchar(GOtermIDs)!=10)){warning('plotGOgraph: Number of characters in GOtermIDs is not correct. Must be 10.')}
if(any(!is.character(GOtermIDs))){stop('plotGOgraph: GOtermIDs required as input. Check that it is IDs and not GOterm numbers. Function stops.')}
# PlotFile Text?
if(!is.character(PlotFile)){stop('plotGOgraph: PlotFile is not of type character! Function stops.')}
# PlotDirectory Text?
if(!is.character(PlotDirectory)){stop('plotGOgraph: PlotDirectory is not of type character! Function stops.')}
# Adj keine Null-Matrix?
if(nrow(Adj)==0){print('plotGOgraph: Adjacency matrix is <0x0 matrix> no data to draw.')
								return()}


FileParts <- fileparts(PlotFile)
# Herausfinden, welches Format die Ausgabedatei haben soll.
if(length(FileParts$name)==0){
	FileName <- 'plotGOgraphOutput'
}else{
	FileName <- FileParts$name
}#end if(length(FileParts$name)==0)

if(nchar(FileParts$ext)==0){
	Ext <- '.png'
}else{
	Ext <- FileParts$ext
}#end if(nchar(FileParts$ext)==0)

MAXLENSTR = 40; # maximale Laenge der GOtermNames in den Knoten

# Parameter zurechtpfriemeln
if(!is.null(GOtermString)){
	ZuLang <- which(nchar(GOtermString)>40)
	Names <- strtrim(GOtermString,MAXLENSTR) # strings auf MAXLENSTR  kuerzen
	Names[ZuLang] <- paste0(strtrim(Names, MAXLENSTR-3), '...')[ZuLang]
	if(nrow(Adj)!=length(GOtermString)){stop('plotGOgraph: Length of "GOtermString" differs from row and column dimension of adjacency matrix. Function stops.')}
}else{
	#require(GO.db)
	#requireNamespace(package ='GO.db', quietly = TRUE)
	GOtermString <- termDescription(GOtermIDs)
	ZuLang <- which(nchar(GOtermString)>40)
	Names <- strtrim(GOtermString,MAXLENSTR) # strings auf MAXLENSTR  kuerzen
	Names[ZuLang] <- paste0(strtrim(Names, MAXLENSTR-3), '...')[ZuLang]
}#end if(!is.null(GOtermString))

DotFile  <- paste0(FileName, '.dot');
RemGiven <- !is.null(Remarkable)
PvaluesGiven = !is.null(Pvalues) ; # Pvalues gegeben?
NrGenesGiven = !is.null(NrGenesInTerm) ;
ExpGiven = !is.null(Expected)
ObsGiven = !is.null(Observed)
ImpGiven = !is.null(Importance)
UpGiven = !is.null(Up); 
AdjPlusEltern <- Adj
if(is.null(Significant)){
	Significant <- rep(0,nrow(Adj))
}
# Significant genauso lang wie nrow(Adj)?
if(nrow(Adj)!=length(Significant)){stop('plotGOgraph: Length of "Significant" differs from row and column dimension of adjacency matrix. Function stops.')}
if(is.null(IsHeadline)){
	IsHeadline <- rep(0,nrow(Adj))
}
# IsHeadline genauso lang wie nrow(Adj)?
if(nrow(Adj)!=length(IsHeadline)){stop('plotGOgraph: Length of "IsHeadline" differs from dimension of adjacency matrix. Function stops.')}


# if(RemGiven){
	# remarkable <- paste0('Rem = ',signif(Remarkable,digits=2));   #remarkable gerundet
	# if(nrow(Adj)!=length(Remarkable)){stop('plotGOgraph: Length of "Remarkable" differs from dimension of adjacency matrix. Function stops.')}
# }#end if (RemGiven)

# ########################################################
# Zusammenbau des Textes, der in die Knoten soll
IndKomplettText <- which(Significant != 0)

Text<- paste0('"',GOtermIDs,'" [label="',GOtermIDs ,'\\n',Names) #dieser Erste Teil ist bei allen Knoten gleich.
#jetzt noch bedingt Text in Knoten einfuegen:
#Teste ob Anzahl Gene in Term gegeben sind, wenn ja, schreibe mit in den Text, der im Knoten stehen soll.
if(NrGenesGiven){ 
	# Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n n = ', NrGenesInTerm[IndKomplettText]);
	Text[which(!is.nan(NrGenesInTerm))] <- paste0(Text[which(!is.nan(NrGenesInTerm))],'\\n n = ', NrGenesInTerm[which(!is.nan(NrGenesInTerm))])
	# Text <- paste0(Text,'\\n n = ', NrGenesInTerm)
}# end if NrGenesGiven
if(ExpGiven){
	EXP <- paste0('Exp = ',signif(Expected,digits=2)) # Expected auf eine NKS runden
	Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n', EXP[IndKomplettText]);
}# end if ExpGiven
if(ObsGiven){ 
	Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n Obs = ', Observed[IndKomplettText]);
}# end if ObsGiven
#Teste ob P-values uebergeben wurden, wenn ja, schreibe mit in den Text, der im Knoten stehen soll.
if(PvaluesGiven){
	PV <- paste0('Pval = ',signif(Pvalues,digits=2)) #runde p-Values auf 1 NKS und schreibe "Pval= " davor
	Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n',PV[IndKomplettText])# Vektor mit Text, der in Knoten stehen soll, "die gegeben" waren (1.Teil Knotentext) mit P-value
}# end if PvaluesGiven
if(ImpGiven){ 
	Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n Imp = ', Importance[IndKomplettText]);
}# end if ImpGiven
if(RemGiven){
Text[IndKomplettText] <- paste0(Text[IndKomplettText],'\\n Rem = ', Remarkable[IndKomplettText])#  Remarkable Wert ganz unten im Knoten
}# end if RemGiven
Text <- paste0(Text,'"];') 

# Text fertig zusammengebaut.
# ########################################################

# ########################################################
# Bedingt anpassen welche Knoten welche Farbe haben sollen:
# Der Rand zeigt immer die Knotenfarbe an (ROT/GRUEN je nach Up).
# Blaetter die nicht Headlines sind: Blau mit Rand in der Knotenfarbe (-> Variablen: BlauRot, BlauGruen)
# Blaetter die Headlines sind: Gelb mit Rand in der Knotenfarbe, Schrift blau (-> Variablen: BlauGelbRot, BlauGelbGruen)
# Innere Knoten die Headlines sind: Gelb mit Rand in der Knotenfarbe (-> Variablen: RotGelb, GruenGelb)
# Signifikante Knoten (ohne sonstige Eigentschaft): Rot/Gruen je nach Up (-> Variablen: Rot, Gruen)
# Nicht-Signifikante Knoten weiss

Rot<-''
RotGelb<-''
Gruen<-''
GruenGelb<-''
Blau <- ''
BlauRot<-''
BlauGruen<-''
BlauGelbRot<-''
BlauGelbGruen<-''
if(length(GOtermIDs[as.logical(Significant)])>0){
	Rot <- paste0('"', GOtermIDs[as.logical(Significant)], '" [color="red", style=filled];')
}#end if signifikante knoten rot
if(MarkDetails & (length(Matrix::colSums(Adj) == 0)>0)){
	Blau <- paste0('"', GOtermIDs[Matrix::colSums(Adj) == 0], '" [color="lightblue2", style=filled];')
}# end if Blaetter (egal ob Significant oder nicht!)
if(length(GOtermIDs[Significant==1&IsHeadline==1])>0){
	RotGelb <- paste0('"',GOtermIDs[Significant==1&IsHeadline==1],'" [color="red", fillcolor="#FFFF80",penwidth=3];')
}#end if: Knoten, die Headlines sind, werden von der Farbe umrandet, die der Knoten zuvor hatte. (hier rot)
if(UpGiven){ #nur die, die signifikant und schon rot sind, koennen gruen werden. und zwar die, wo Up==0 ist.
	if(length(GOtermIDs[Significant==1&Up==0])>0){
		Gruen<-paste0('"',GOtermIDs[Significant==1&Up==0],'" [color="green", style=filled];')
	}#end if: gruene Knoten gefaerbt
	if(length(GOtermIDs[Significant==1&Up==0&IsHeadline==1])>0){
		GruenGelb <- paste0('"',GOtermIDs[Significant==1&Up==0&IsHeadline==1],'" [fillcolor="#FFFF80", color="green",penwidth=3];')
	}#end if: Knoten, die Headlines sind, werden von der Farbe umrandet, die der Knoten zuvor hatte. (hier gruen)
	if(MarkDetails){ #falls also gewuenscht ist, dass die Blaetter blau gefaerbt werden.
		if(length(GOtermIDs[Significant==1&IsHeadline==0&Up==1&(Matrix::colSums(Adj) == 0)])>0){
			BlauRot <- paste0('"',GOtermIDs[IsHeadline==0&Up==1&(Matrix::colSums(Adj) == 0)],'" [fillcolor="lightblue2", color="red", penwidth=3];')
		}# end if: Blaetter Blau mit Rotem Rand
		if(length(GOtermIDs[Significant==1&IsHeadline==0&Up==0&(Matrix::colSums(Adj) == 0)])>0){
			BlauGruen<-paste0('"',GOtermIDs[IsHeadline==0&Up==0&(Matrix::colSums(Adj) == 0)],'" [fillcolor="lightblue2", color="green", penwidth=3];')
		}# end if: Blaetter Blau mit Gruenem Rand
		if(length(GOtermIDs[Significant==1&IsHeadline==1&Up==1&(Matrix::colSums(Adj) == 0)])>0){
			BlauGelbRot<-paste0('"',GOtermIDs[IsHeadline==1&Up==1&(Matrix::colSums(Adj) == 0)],'" [fillcolor="#FFFF80", color="red", fontcolor="blue", penwidth=3];')
		}# end if: Blaetter Blaue Schrift bei Headlines (gelb) und umrandet mit rot 
		if(length(GOtermIDs[Significant==1&IsHeadline==1&Up==0&(Matrix::colSums(Adj) == 0)])>0){
			BlauGelbGruen<-paste0('"',GOtermIDs[IsHeadline==1&Up==0&(Matrix::colSums(Adj) == 0)],'" [fillcolor="#FFFF80", color="green", fontcolor="blue", penwidth=3];')
		}#end if: Blaetter Blaue Schrift bei Headlines (gelb) und umrandet mit gruen
	}#end if: MarkDetails - Blaue Blaetter
}else{ #wenn Up nicht gegeben ist aber MarkDetails, dann mit rotem Rand Blaetter blau, aber auch nur fuer Signifikante Knoten
	if(MarkDetails){ 
		if(length(GOtermIDs[Significant==1&IsHeadline==0&(Matrix::colSums(Adj) == 0)])>0){
			BlauRot<-paste0('"',GOtermIDs[IsHeadline==0&(Matrix::colSums(Adj) == 0)],'" [fillcolor="lightblue2", color="red", penwidth=3];')
			}#end if: Blaetter Blau mit rotem Rand
		if(length(GOtermIDs[Significant==1&IsHeadline==1&(Matrix::colSums(Adj) == 0)])>0){
			BlauGelbRot<-paste0('"',GOtermIDs[IsHeadline==1&(Matrix::colSums(Adj) == 0)],'" [fillcolor="#FFFF80", color="red", fontcolor="blue", penwidth=3];')
			}# end if: Blaetter, die Headlines sind gelb mit blauer Schrift und rotem Rand
	}#end if MarkDetails
}#end if+else UpGiven
# Knoten haben jetzt die richtigen Farben
# ########################################################

# ########################################################
# Jetzt noch die Verbindungen der Knoten untereinander aufschreiben: Vielleicht fuer bessere Laufzeit for-Schleifen eliminieren??
count<-1
Edges<-c('')
for (i in c(1:length(GOtermIDs))){
	for (j in c(1:length(GOtermIDs))){ # Adjazenzmatrix ist quadratisch => nrow(AdjPlusEltern)== ncol(AdjPlusEltern)
		if (AdjPlusEltern[i,j]==1){
			Edges[count]<-paste0('"',GOtermIDs[j],'" -> "',GOtermIDs[i],'"')
			count<-count+1
		}#end if
	}#end for (j)
}# end for (i)
# alle Verbindungen sind in Dot-Format in 'Edges' gespeichert
# ######################################################### 

# #########################################################
#jetzt dot-Code in eine Datei ausschreiben, Konsole aufrufen dot-file ausfuehren und png erstellen, dot-File loeschen

CurrentDir = getwd(); # merken in welchem Verzeichnis wir grade sind
setwd(PlotDirectory)
dot <- c('digraph G{',Text,Rot,RotGelb,Gruen,GruenGelb,Blau,BlauRot,BlauGruen,BlauGelbRot,BlauGelbGruen,Edges,'}') 
write(dot, file = DotFile) #dot-File erstellen

# Schaue, ob Output-Filename schon vergeben ist, wenn ja, fuege fortlaufende Nummer hinten an:
if(!Overwrite){
	VV = checkFilename(Directory = PlotDirectory, FileName = FileName, Extension = Ext, ReadOrWrite = FALSE, NameOfFunctionCalled = 'plotGOgraph')
}else{
	VV <- list(Filename = paste0(FileName,Ext))
}# end if(!Overwrite)

oeffnen <- paste0("dot -T", gsub('\\.','',Ext)," ", DotFile," -o ", VV$Filename," -Gsize=900,900 -Nfontname=arial") #muesste auf Linux und Windows laufen
#Hinweis: wenn man in oeffnen noch "-Nfontname=helvetica ..." einfuegt, kann man noch die Schriftart aendern (... ersetzen ohne "")
try(system(oeffnen)) #dot-File oeffnen und ausfuehren
OSname <- Sys.info()['sysname']
# *.dot file loeschen je nach Betriebssystem
switch(OSname,  #checke ob Linux oder Windows oder sonst ein Betriebssystem
	Windows = {	loeschen <- paste0("del ",DotFile)
							try(shell(loeschen))#dot-File gleich wieder loeschen. Es wird nur das (png-)File gespeichert.
						},
	Linux = {	loeschen <- paste0("rm ",DotFile)
						try(system(loeschen))#dot-File gleich wieder loeschen. Es wird nur das (png-)File gespeichert.
					},
	Darwin = 	{	print("Mac OS: *.dot file has to be removed by hand.")
						},
	'Unknown OS. *.dot file has to be removed by hand.'
)#end switch OSname

print(paste0('plotGOgraph: ',gsub('\\.','',Ext),'-File named "',VV$Filename,'" saved in ',PlotDirectory)); # gib laut 
setwd(CurrentDir);     # wieder auf CurrentDir zuruecksetzen

} # end function plotGOgraph
 