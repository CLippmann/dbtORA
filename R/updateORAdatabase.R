updateORAdatabase <- function(){
# Function to update the ORA database and save in ORA/inst/extdata

# OUTPUT:
# Following files will be updated:
# AllAnnNCBIsPlusGeneName.names			Contains all NCBIs that are annotated to at least one GOterm.
# AdjBPsparseMatrix.lrn							Sparse Matrix of BP GOterms. Represents DAG structure.
# AdjBPsparseMatrix.names						Names of BP GOterms corresponding to AdjBPsparseMatrix.
# AdjMFsparseMatrix.lrn							Sparse Matrix of MF GOterms. Represents DAG structure.
# AdjMfsparseMatrix.names						Names of MF GOterms corresponding to AdjMFsparseMatrix.
# AdjCCsparseMatrix.lrn							Sparse Matrix of CC GOterms. Represents DAG structure.
# AdjCCsparseMatrix.names						Names of CC GOterms corresponding to AdjCCsparseMatrix.
# GOAall.lrn												NCBIs and GOtermNrs they are annotated to as adjacency list.
# GOTermInfosBP.lrn									Statistical infos about BP GOterms.
# GOTermInfosMF.lrn									Statistical infos about MF GOterms.	
# GOTermInfosCC.lrn									Statistical infos about CC GOterms.
# BPAncestorsSparseMatrix.lrn				Sparse Matrix of BP GOterms and all their ancestors.
# BPAncestorsSparseMatrix.names			Names of BP GOterms corresponding to BPAncestorsSparseMatrix.
# MFAncestorsSparseMatrix.lrn				Sparse Matrix of MF GOterms and all their ancestors.
# MFAncestorsSparseMatrix.names			Names of MF GOterms corresponding to MFAncestorsSparseMatrix.
# CCAncestorsSparseMatrix.lrn				Sparse Matrix of CC GOterms and all their ancestors.
# CCAncestorsSparseMatrix.names			Names of CC GOterms corresponding to CCAncestorsSparseMatrix.

print('Loading packages.')

##### Lade wichtigste Datenbanken: ####
print('Loading data.')
GO2EG <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egGO2EG) # Fuer direkte Anns
GO2ALLEGS <- AnnotationDbi::as.list(org.Hs.eg.db::org.Hs.egGO2ALLEGS) # Fuer direkte und indirekte Anns
BPchildren <- AnnotationDbi::as.list(GO.db::GOBPCHILDREN) # Kinder in BP fuer DAG-Struktur
MFchildren <- AnnotationDbi::as.list(GO.db::GOMFCHILDREN) # Kinder in MF fuer DAG-Struktur
CCchildren <- AnnotationDbi::as.list(GO.db::GOCCCHILDREN) # Kinder in CC fuer DAG-Struktur
GOBPAnc <- AnnotationDbi::as.list(GO.db::GOBPANCESTOR) # Vorfahren in BP fuer DAG-Struktur
GOMFAnc <- AnnotationDbi::as.list(GO.db::GOMFANCESTOR) # Vorfahren in MF fuer DAG-Struktur
GOCCAnc <- AnnotationDbi::as.list(GO.db::GOCCANCESTOR) # Vorfahren in CC fuer DAG-Struktur

##### Allgemeine Parameter: #####
OutDirectory <- system.file('extdata', package = 'ORA')

##########################################################################################
## GOAall.lrn ############################################################################
##########################################################################################

# Aus der Dokumentation:
# org.Hs.egGO is an R object that provides mappings between entrez gene identifiers and the GO
# identifiers that they are directly associated with.  This mapping and its reverse mapping do NOT
# associate the child terms from the GO ontology with the gene.  Only the directly evidenced terms
# are represented here.
# org.Hs.egGO2ALLEGS is an R object that provides mappings between a given GO identifier and
# all of the Entrez Gene identifiers annotated at that GO term OR TO ONE OF ITS CHILD NODES
# in the GO ontology. Thus, this mapping is much larger and more inclusive than org.Hs.egGO2EG

# als Liste:
AnzNCBIProGOTerm <- sapply(GO2ALLEGS,length)
NCBIs <- unlist(unname(GO2ALLEGS))
GOTerme <- rep(names(AnzNCBIProGOTerm),AnzNCBIProGOTerm)
Evidence <- names(NCBIs)
OntologyName <- Ontology(GOTerme)
OntoNameNAind <- which(is.na(OntologyName))
if(length(OntoNameNAind)!=0){
NCBIs <- NCBIs[-OntoNameNAind]
GOTerme <- GOTerme[-OntoNameNAind]
Evidence <- Evidence[-OntoNameNAind]
OntologyName <- OntologyName[-OntoNameNAind]
}
Key <- seq(1,length(NCBIs))
FileName <- paste0('GOAall.lrn')
TempMatrix <- cbind(as.numeric(NCBIs), termNr(GOTerme), EvidenceCode2Key(Evidence), OntoName2No(OntologyName))
temp2 <- sortby(TempMatrix, 2)
Data <- sortby(temp2, 1)
Header <- c('Key',	'AnnNCBI',	'GOTermNo',	'EvidenceKey',	'OntologyNo')
Comments <- paste0('This file contains infos from "org.Hs.egGO2ALLEGS". \n# From the documentation of the package:\n# org.Hs.egGO2ALLEGS is an R object that provides mappings between a given GO identifier and\n# all of the Entrez Gene identifiers annotated at that GO term OR TO ONE OF ITS CHILD NODES\n# in the GO ontology. Thus, this mapping is much larger and more inclusive than org.Hs.egGO2EG.\n# Last update: ', Sys.Date(),'.')
# WriteLRN(FileName,Data,Header,Key,DataDefined,OutDirectory,CommentOrDigits);
suppressWarnings(WriteLRN(FileName,Data,Header,Key,c(),OutDirectory,Comments));
print('Generate updatet files.')
print('1/18 done.')


##########################################################################################
##########################################################################################
# Schreibe Struktur der GeneOntology aus GO.db aus.
# In Form einer sparse Matrix.

##########################################################################################
## AdjCCsparseMatrix.lrn und ~.names######################################################
##########################################################################################
# Zuerst sortieren wir mal die Liste nach ihren Namen (also den GO-Termen):
n <- termNr(names(CCchildren)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
CClistSort <- CCchildren[o]
# leere Zeilen:
IstBlattCC <- is.na(CClistSort)
# Anz Kinder pro Term
AnzKinderCCInNichtNullZeilen <- sapply(CClistSort[!IstBlattCC], length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Kinder in den Nicht-Null-Zeilen gibt.
ZeilenCC <- rep(which(!IstBlattCC), AnzKinderCCInNichtNullZeilen)
SpaltenCC <- match(unname(unlist(unname(CClistSort[!IstBlattCC]))),termId(orderedTerms))
InhaltCC <- 1
AdjMatrixCC <- sparseMatrix(ZeilenCC,SpaltenCC,x=InhaltCC, dims=c(length(orderedTerms), length(orderedTerms)), index1=TRUE, dimnames=list(orderedTerms,orderedTerms))
Comment <- 'Saving sparseMatrix of all GO-Terms in CC. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] is parent of GOTerms[SpaltenInd[i]]. The corresponding \n# GOTerms can be found in AdjCCsparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('AdjCCsparseMatrix', ZeilenCC, SpaltenCC, InhaltCC, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('3/18 done.')

##########################################################################################
## AdjMFsparseMatrix.lrn und ~.names######################################################
##########################################################################################
# Zuerst sortieren wir mal die Liste nach ihren Namen (also den GO-Termen):
n <-termNr(names(MFchildren)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
MFlistSort <- MFchildren[o]
# leere Zeilen:
IstBlattMF <- is.na(MFlistSort)
# Anz Kinder pro Term
AnzKinderMFInNichtNullZeilen <- sapply(MFlistSort[!IstBlattMF], length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Kinder in den Nicht-Null-Zeilen gibt.
ZeilenMF <- rep(which(!IstBlattMF), AnzKinderMFInNichtNullZeilen)
SpaltenMF <- match(unname(unlist(unname(MFlistSort[!IstBlattMF]))),termId(orderedTerms))
InhaltMF <- 1
AdjMatrixMF <- sparseMatrix(ZeilenMF,SpaltenMF,x=InhaltMF, dims=c(length(orderedTerms), length(orderedTerms)), index1=TRUE, dimnames=list(orderedTerms,orderedTerms))
Comment <- 'Saving sparseMatrix of all GO-Terms in MF. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] is parent of GOTerms[SpaltenInd[i]]. The corresponding \n# GOTerms can be found in AdjMFsparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('AdjMFsparseMatrix', ZeilenMF, SpaltenMF, InhaltMF, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('5/18 done.')

##########################################################################################
## AdjBPsparseMatrix.lrn und ~.names######################################################
##########################################################################################
# Zuerst sortieren wir mal die Liste nach ihren Namen (also den GO-Termen):
n <-termNr(names(BPchildren)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
BPlistSort <- BPchildren[o]
# leere Zeilen:
IstBlattBP <- is.na(BPlistSort)
# Anz Kinder pro Term
AnzKinderBPInNichtNullZeilen <- sapply(BPlistSort[!IstBlattBP], length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Kinder in den Nicht-Null-Zeilen gibt.
ZeilenBP <- rep(which(!IstBlattBP), AnzKinderBPInNichtNullZeilen)
SpaltenBP <- match(unname(unlist(unname(BPlistSort[!IstBlattBP]))),termId(orderedTerms))
InhaltBP <- 1
AdjMatrixBP <- sparseMatrix(ZeilenBP,SpaltenBP,x=InhaltBP, dims=c(length(orderedTerms), length(orderedTerms)), index1=TRUE, dimnames=list(orderedTerms,orderedTerms))
Comment <- 'Saving sparseMatrix of all GO-Terms in BP. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] is parent of GOTerms[SpaltenInd[i]]. The corresponding \n# GOTerms can be found in AdjBPsparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('AdjBPsparseMatrix', ZeilenBP, SpaltenBP, InhaltBP, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('7/18 done.')





##########################################################################################
##########################################################################################
# zum Speichern von Informationen aus GOBPANCESTOR, GOMFANCESTOR und GOCCANCESTOR aus dem Package GO.db.

##########################################################################################
## BPAncestorsSparseMatrix.lrn und ~.names ###############################################
##########################################################################################
GOBPAncWOall <- lapply(GOBPAnc, function(GOi) GOi[GOi!='all']) # "all" als Vorfahre entfernen
n <-termNr(names(GOBPAncWOall)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
GOBPAncSort <- GOBPAncWOall[o]
# Anz Vorfahren pro Term
AnzVorfahren <- sapply(GOBPAncSort, length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Vorfahren gibt.
ZeilenBP <- match(rep(names(GOBPAncSort), AnzVorfahren), names(GOBPAncSort))
SpaltenBP <- match(unlist(unname(GOBPAncSort)), names(GOBPAncSort))
InhaltBP <- 1
#AncAdjMatrixBP <- sparseMatrix(ZeilenBP, SpaltenBP, x = InhaltBP, dims = c(length(orderedTerms), length(orderedTerms)), index1 = TRUE, dimnames = list(termId(orderedTerms), termId(orderedTerms)))
Comment <- 'Saving sparseMatrix of all GO-Terms and their ancestors in BP. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] has GOTerms[SpaltenInd[i]] as ancestor. The corresponding \n# GOTerms can be found in BPAncestorsSparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('BPAncestorsSparseMatrix', ZeilenBP, SpaltenBP, InhaltBP, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('9/18 done.')

##########################################################################################
## MFAncestorsSparseMatrix.lrn und ~.names ###############################################
##########################################################################################
GOMFAncWOall <- lapply(GOMFAnc, function(GOi) GOi[GOi!='all']) # "all" als Vorfahre entfernen
n <-termNr(names(GOMFAncWOall)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
GOMFAncSort <- GOMFAncWOall[o]
# Anz Vorfahren pro Term
AnzVorfahren <- sapply(GOMFAncSort, length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Vorfahren gibt.
ZeilenMF <- match(rep(names(GOMFAncSort), AnzVorfahren), names(GOMFAncSort))
SpaltenMF <- match(unlist(unname(GOMFAncSort)), names(GOMFAncSort))
InhaltMF <- 1
#AncAdjMatrixMF <- sparseMatrix(ZeilenMF, SpaltenMF, x = InhaltMF, dims = c(length(orderedTerms), length(orderedTerms)), index1 = TRUE, dimnames = list(termId(orderedTerms), termId(orderedTerms)))
Comment <- 'Saving sparseMatrix of all GO-Terms and their ancestors in MF. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] has GOTerms[SpaltenInd[i]] as ancestor. The corresponding \n# GOTerms can be found in MFAncestorsSparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('MFAncestorsSparseMatrix', ZeilenMF, SpaltenMF, InhaltMF, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('11/18 done.')

##########################################################################################
## CCAncestorsSparseMatrix.lrn und ~.names ###############################################
##########################################################################################
GOCCAncWOall <- lapply(GOCCAnc, function(GOi) GOi[GOi!='all']) # "all" als Vorfahre entfernen
n <-termNr(names(GOCCAncWOall)) # Wir nehmen termNo, damit wir 1,2,3.. und nicht 1,10,111,... sortieren
o <- order(n)
orderedTerms <- n[o]
GOCCAncSort <- GOCCAncWOall[o]
# Anz Vorfahren pro Term
AnzVorfahren <- sapply(GOCCAncSort, length)
# Fuer SparseMatrix brauchen wir die Indizes der nicht-Null-Eintraege.
# Nicht-Null-Zeilen so oft wiederholt, wie es Vorfahren gibt.
ZeilenCC <- match(rep(names(GOCCAncSort), AnzVorfahren), names(GOCCAncSort))
SpaltenCC <- match(unlist(unname(GOCCAncSort)), names(GOCCAncSort))
InhaltCC <- 1
#AncAdjMatrixCC <- sparseMatrix(ZeilenCC, SpaltenCC, x = InhaltCC, dims = c(length(orderedTerms), length(orderedTerms)), index1 = TRUE, dimnames = list(termId(orderedTerms), termId(orderedTerms)))
Comment <- 'Saving sparseMatrix of all GO-Terms and their ancestors in CC. Each row describes a connection between Terms. \n# Meaning GOTerms[ZeilenInd[i]] has GOTerms[SpaltenInd[i]] as ancestor. The corresponding \n# GOTerms can be found in CCAncestorsSparseMatrix.names.'
# WriteSparseMatrix(FileName, ZeilenInd, SpaltenInd, Inhalt, Dimension, DimNames, OutDirectory, Header, Key, CommentOrDigits){
suppressWarnings(WriteSparseMatrix('CCAncestorsSparseMatrix', ZeilenCC, SpaltenCC, InhaltCC, c(length(orderedTerms),length(orderedTerms)), list(termId(orderedTerms),termId(orderedTerms)), OutDirectory, c('Key', 'ZeilenInd', 'SpaltenInd', 'Inhalt'), Comment=Comment))
print('13/18 done.')





##########################################################################################
##########################################################################################
# Pro Ontologie eine *.lrn mit den wichtigsten Daten zu den 
# GO-Termen in dieser Ontologie. Dazu gehoeren jeweils:

# LRN:
# Pro GO-Term:
# Key															= GO-Term Nummer
# AnzGeneDirektAnnotiert 					= Anzahl der Gene, die direkt an diesem Knoten annotiert sind.
# AlleAnnotierungen 							= Anzahl der Gene, die direkt oder indirekt, d.h. an Nachkommen dieses Knotens, annotiert sind.
# AnzGeneDirektAnnotiertMan 			= Anzahl der Gene, die direkt an diesem Knoten annotiert sind ohne "IEA". Also nur manually curated.
# AnzGeneDirektAnnotiertAllUniq 	= Anzahl der Gene, die direkt an diesem Knoten annotiert sind ohne den Evidence-Code zu beachten. Also unique.
# AnzGeneDirektAnnotiertManUniq 	= Anzahl der Gene, die direkt an diesem Knoten annotiert sind und manually curated. Nach der Auswahl der manually curated Annotierungen wird ein Gen mit unterschiedlichen Evidence-Codes nur einmal gezaehlt. Also unique.
# AlleAnnotierungenManuUniq				= Anzahl der Gene, die direkt oder indirekt annotiert sind ohne IEA also nur manually curated.
# AnzahlNachkommen 								= Anzahl der Nachkommen dieses Terms, also Kinder und Kindeskinder usw..
# AnzahlKinder 										= Anzahl der Kinder des Terms, d.h. nur direkte Nachkommen.
# AnzahlVorfahren 								= Anzahl der Vorfahren des Terms bis zur Wurzel.
# AnzahlEltern 										= Anzahl der Eltern, d.h. der direkten Vorfahren.
# InformationContent 							= negativer Logarithmus des Verhaeltnisses von der Anzahl der unique, indirekten und direkten Annotierungen des Terms zu allen unique, indirekten und direkten Annotierungen in der Ontologie (letzteres ist die Anzahl der indirekten und direkten Annotierungen der Wurzel).
# Tiefe 													= Laenge des kuerzesten Pfads vom Knoten zur Wurzel.
# Level 													= Laenge des laengsten Pfads vom Knoten zur Wurzel.
# Hoehe 													= Laenge des kuerzesten Pfads vom Knoten zu einem seiner Blaetter.
# HoehenLevel 										= Laenge des laengsten Pfads vom Knoten zu einem seiner Blaetter.

# erstelle Kommentar fuer AusgabeDateien
Comments = paste('Infos about GOterms. Last update: ', Sys.Date())
AnzGeneDirektAnnotiert <- sapply(GO2EG,length) 
GeneDirektAnnotiertMan <- lapply(GO2EG, function(GO2EGi){if(any(names(GO2EGi)=='IEA')){GO2EGi <- GO2EGi[-which(names(GO2EGi) == 'IEA')]}else{return(GO2EGi)}}) #nur ManuCur Gene
AnzGeneDirektAnnotiertMan <- sapply(GeneDirektAnnotiertMan,length) 
AnzGeneDirektAnnotiertAllUniq <- sapply(lapply(GO2EG,unique),length) #nur unique 
GeneDirektAnnotiertManUniq <- lapply(GeneDirektAnnotiertMan,unique) # Manu und unique Gene
AnzGeneDirektAnnotiertManUniq <- sapply(GeneDirektAnnotiertManUniq,length)
#Indirekte Ann. fuer Information Content.
AnzAllAnnotations <- sapply(lapply(GO2ALLEGS, unique), length)
AnzManuAnnotations <- sapply(lapply(lapply(GO2ALLEGS, function(GO2EGi){if(any(names(GO2EGi)=='IEA')){GO2EGi <- GO2EGi[-which(names(GO2EGi) == 'IEA')]}else{return(GO2EGi)}}), unique), length)

##########################################################################################
## GOTermInfosBP.lrn #####################################################################
##########################################################################################
zz <- as.list(GO.db::GOBPOFFSPRING)
AnzahlNachkommenBP <- sapply(zz,length)
AnzahlNachkommenBP[which(is.na(zz))] <- 0
AnzahlKinderBP <- sapply(BPchildren,length)
AnzahlKinderBP[which(is.na(BPchildren))] <- 0
AnzahlVorfahrenBP <- sapply(GOBPAnc,length)
AnzahlVorfahrenBP[which(GOBPAnc=='all')] <- 0
yy <- as.list(GO.db::GOBPPARENTS)
KeyBP <- termNr(names(yy)) # Key = GO-Term Nummern aus yy, da in GO2EG nur diejenigen sind, die mindestens eine direkte Annotierung haben. Hier nur BP. (yy ist viel groesser)
AnzahlElternBP <- sapply(yy,length)
AnzahlElternBP[which(yy=='all')] <- 0
IndexElternBP <- match(names(GO2EG),names(yy))
NichtNA <- which(!is.na(IndexElternBP)) #Welche von GO2EG sind in yy
IndexElternBP <- IndexElternBP[NichtNA]
AnzBPTerme <- length(KeyBP)
AnzGeneDirektAnnotiertBP <- c(rep(0, AnzBPTerme))
AnzGeneDirektAnnotiertManBP <- c(rep(0, AnzBPTerme))
AnzGeneDirektAnnotiertAllUniqBP <- c(rep(0, AnzBPTerme))
AnzGeneDirektAnnotiertManUniqBP <- c(rep(0, AnzBPTerme))
AnzGeneDirektAnnotiertBP[IndexElternBP] <- AnzGeneDirektAnnotiert[NichtNA] 
AnzGeneDirektAnnotiertManBP[IndexElternBP] <- AnzGeneDirektAnnotiertMan[NichtNA] 
AnzGeneDirektAnnotiertAllUniqBP[IndexElternBP] <- AnzGeneDirektAnnotiertAllUniq[NichtNA]
AnzGeneDirektAnnotiertManUniqBP[IndexElternBP] <- AnzGeneDirektAnnotiertManUniq[NichtNA]
BPAnnotInd <- match(termId(KeyBP),names(AnzAllAnnotations))
AnzAllAnnotationsBP <- KeyBP*0
AnzAllAnnotationsBP[which(!is.na(BPAnnotInd))] <- AnzAllAnnotations[BPAnnotInd[which(!is.na(BPAnnotInd))]]
Pti <-  AnzAllAnnotationsBP/AnzAllAnnotations[which(names(AnzAllAnnotations)=='GO:0008150')] # Annot. an Term unique/ alle annotierten Gene an Wurzel
BPAnnotManuInd <- match(termId(KeyBP), names(AnzManuAnnotations))
AnzAllManuAnnotationsBP <- KeyBP*0
AnzAllManuAnnotationsBP[which(!is.na(BPAnnotManuInd))] <- AnzManuAnnotations[BPAnnotManuInd[which(!is.na(BPAnnotManuInd))]]
PtiManu <- AnzAllManuAnnotationsBP/AnzManuAnnotations[which(names(AnzManuAnnotations)=='GO:0008150')]
#Minimum ist 1/16655 = 6.004203e-05
#Dementsprechend Maximum InformationContent = -log2(1/16655) = 14.02367. (gerundet 14.02)
#Idee: Evtl sollte man das zur Vergleichbarkeit mit anderen Ontologien in Prozent von Maximum angeben... 
NullIndex <- which(Pti==0)
InformationContentBP <- round(x = -log2(Pti), digits = 2)
InformationContentBP[NullIndex] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = NaN.
NullIndexManu <- which(PtiManu==0)
InformationContentBPManu <- round(x = -log2(PtiManu), digits = 2)
InformationContentBPManu[NullIndexManu] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = NaN.
AlleAnnotierungenBP <- AnzAllAnnotationsBP
TiefeBP <- MinPathL2Root('BP')
LevelBP <- MaxPathL2Root('BP')
HoeheBP <- MinPathL2Leaf('BP')
HoehenLevelBP <- MaxPathL2Leaf('BP')

# erstelle aus Daten Matrix
Data <- cbind(KeyBP, AnzGeneDirektAnnotiertBP, AlleAnnotierungenBP, AnzGeneDirektAnnotiertManBP, AnzGeneDirektAnnotiertAllUniqBP, AnzGeneDirektAnnotiertManUniqBP, AnzAllManuAnnotationsBP, AnzahlNachkommenBP, AnzahlKinderBP, AnzahlVorfahrenBP, AnzahlElternBP, InformationContentBP, InformationContentBPManu, TiefeBP, LevelBP, HoeheBP, HoehenLevelBP)
#Schreibe aus.
#WriteLRN(FileName = ,Data = ,Header = ,Key = ,DataDefined = ,OutDirectory = ,CommentOrDigits =)
suppressWarnings(WriteLRN('GOTermInfosBP', Data = Data, Header = c('GOtermNr', 'AnzGeneDirektAnnotiertBP', 'AlleAnnotierungenUniqBP', 'AnzGeneDirektAnnotiertManBP', 'AnzGeneDirektAnnotiertAllUniqBP', 'AnzGeneDirektAnnotiertManUniqBP', 'AnzAllManuAnnotationsBP', 'AnzahlNachkommenBP', 'AnzahlKinderBP', 'AnzahlVorfahrenBP', 'AnzahlElternBP', 'InformationContentBP', 'InformationContentBPManu', 'MinPathL2RootBP', 'MaxPathL2RootBP', 'MinPathL2LeafBP', 'MaxPathL2LeafBP'), 
Key = KeyBP, DataDefined = c(), OutDirectory = OutDirectory, CommentOrDigits = Comments))
print('14/18 done.')

##########################################################################################
## GOTermInfosMF.lrn #####################################################################
##########################################################################################
zz <- as.list(GO.db::GOMFOFFSPRING)
AnzahlNachkommenMF <- sapply(zz,length)
AnzahlNachkommenMF[which(is.na(zz))] <- 0
AnzahlKinderMF <- sapply(MFchildren,length)
AnzahlKinderMF[which(is.na(MFchildren))] <- 0
AnzahlVorfahrenMF <- sapply(GOMFAnc,length)
AnzahlVorfahrenMF[which(GOMFAnc=='all')] <- 0
yy <- as.list(GO.db::GOMFPARENTS)
KeyMF <- termNr(names(yy)) # Key = GO-Term Nummern aus yy, da in GO2EG nur diejenigen sind, die mindestens eine direkte Annotierung haben. Hier nur MF.
AnzahlElternMF <- sapply(yy,length)
AnzahlElternMF[which(yy=='all')] <- 0
IndexElternMF <- match(names(GO2EG),names(yy))
NichtNA <- which(!is.na(IndexElternMF)) #Welche von GO2EG sind in yy
IndexElternMF <- IndexElternMF[NichtNA]
AnzMFTerme <- length(yy)
AnzGeneDirektAnnotiertMF <- c(rep(0, AnzMFTerme))
AnzGeneDirektAnnotiertManMF <- c(rep(0, AnzMFTerme))
AnzGeneDirektAnnotiertAllUniqMF <- c(rep(0, AnzMFTerme))
AnzGeneDirektAnnotiertManUniqMF <- c(rep(0, AnzMFTerme))
AnzGeneDirektAnnotiertMF[IndexElternMF] <- AnzGeneDirektAnnotiert[NichtNA] 
AnzGeneDirektAnnotiertManMF[IndexElternMF] <- AnzGeneDirektAnnotiertMan[NichtNA] 
AnzGeneDirektAnnotiertAllUniqMF[IndexElternMF] <- AnzGeneDirektAnnotiertAllUniq[NichtNA]
AnzGeneDirektAnnotiertManUniqMF[IndexElternMF] <- AnzGeneDirektAnnotiertManUniq[NichtNA]
MFAnnotInd <- match(termId(KeyMF),names(AnzAllAnnotations))
AnzAllAnnotationsMF <- KeyMF*0
AnzAllAnnotationsMF[which(!is.na(MFAnnotInd))] <- AnzAllAnnotations[MFAnnotInd[which(!is.na(MFAnnotInd))]]
Pti <-  AnzAllAnnotationsMF/AnzAllAnnotations[which(names(AnzAllAnnotations)=='GO:0003674')] # Annot. an Term / Annot. an Wurzel
MFAnnotManuInd <- match(termId(KeyMF), names(AnzManuAnnotations))
AnzAllManuAnnotationsMF <- KeyMF*0
AnzAllManuAnnotationsMF[which(!is.na(MFAnnotManuInd))] <- AnzManuAnnotations[MFAnnotManuInd[which(!is.na(MFAnnotManuInd))]]
PtiManu <- AnzAllManuAnnotationsMF/AnzManuAnnotations[which(names(AnzManuAnnotations)=='GO:0003674')]
#Minimum ist 1/16742 = 5.973002e-05
#Dementsprechend Maximum InformationContent = -log2(1/16742) = 14.03118. (gerundet 14.03)
#Idee: Evtl sollte man das zur Vergleichbarkeit mit anderen Ontologien in Prozent von Maximum angeben... 
NullIndex <- which(Pti==0)
InformationContentMF <- round(x = -log2(Pti), digits = 2)
InformationContentMF[NullIndex] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = -1.
NullIndexManu <- which(PtiManu==0)
InformationContentMFManu <- round(x = -log2(PtiManu), digits = 2)
InformationContentMFManu[NullIndexManu] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = NaN.
AlleAnnotierungenMF <- AnzAllAnnotationsMF
TiefeMF <- MinPathL2Root('MF')
LevelMF <- MaxPathL2Root('MF')
HoeheMF <- MinPathL2Leaf('MF')
HoehenLevelMF <- MaxPathL2Leaf('MF')

# erstelle aus Daten Matrix
Data <- cbind(KeyMF, AnzGeneDirektAnnotiertMF, AnzGeneDirektAnnotiertManMF, AlleAnnotierungenMF, AnzGeneDirektAnnotiertAllUniqMF, AnzGeneDirektAnnotiertManUniqMF, AnzAllManuAnnotationsMF, AnzahlNachkommenMF, AnzahlKinderMF, AnzahlVorfahrenMF, AnzahlElternMF, InformationContentMF, InformationContentMFManu, TiefeMF, LevelMF, HoeheMF, HoehenLevelMF)
#Schreibe aus.
suppressWarnings(WriteLRN('GOTermInfosMF', Data = Data, Header = c('GOtermNr', 'AnzGeneDirektAnnotiertMF', 'AlleAnnotierungenUniqMF', 'AnzGeneDirektAnnotiertManMF', 'AnzGeneDirektAnnotiertAllUniqMF', 'AnzGeneDirektAnnotiertManUniqMF', 'AnzAllManuAnnotationsMF', 'AnzahlNachkommenMF', 'AnzahlKinderMF', 'AnzahlVorfahrenMF', 'AnzahlElternMF', 'InformationContentMF', 'InformationContentMFManu', 'MinPathL2RootMF', 'MaxPathL2RootMF', 'MinPathL2LeafMF', 'MaxPathL2LeafMF'), 
Key = KeyMF, DataDefined = c(), OutDirectory = OutDirectory, CommentOrDigits = Comments))
print('15/18 done.')

##########################################################################################
## GOTermInfosCC.lrn #####################################################################
##########################################################################################
zz <- as.list(GO.db::GOCCOFFSPRING)
AnzahlNachkommenCC <- sapply(zz,length)
AnzahlNachkommenCC[which(is.na(zz))] <- 0
AnzahlKinderCC <- sapply(CCchildren,length)
AnzahlKinderCC[which(is.na(CCchildren))] <- 0
AnzahlVorfahrenCC <- sapply(GOCCAnc,length)
AnzahlVorfahrenCC[which(GOCCAnc=='all')] <- 0
yy <- as.list(GO.db::GOCCPARENTS)
KeyCC <- termNr(names(yy)) # Key = GO-Term Nummern aus yy, da in GO2EG nur diejenigen sind, die mindestens eine direkte Annotierung haben. Hier nur CC.
AnzahlElternCC <- sapply(yy,length)
AnzahlElternCC[which(yy=='all')] <- 0
IndexElternCC <- match(names(GO2EG),names(yy))
NichtNA <- which(!is.na(IndexElternCC)) #Welche von GO2EG sind in yy
IndexElternCC <- IndexElternCC[NichtNA]
AnzCCTerme <- length(yy)
AnzGeneDirektAnnotiertCC <- c(rep(0, AnzCCTerme))
AnzGeneDirektAnnotiertManCC <- c(rep(0, AnzCCTerme))
AnzGeneDirektAnnotiertAllUniqCC <- c(rep(0, AnzCCTerme))
AnzGeneDirektAnnotiertManUniqCC <- c(rep(0, AnzCCTerme))
AnzGeneDirektAnnotiertCC[IndexElternCC] <- AnzGeneDirektAnnotiert[NichtNA] 
AnzGeneDirektAnnotiertManCC[IndexElternCC] <- AnzGeneDirektAnnotiertMan[NichtNA] 
AnzGeneDirektAnnotiertAllUniqCC[IndexElternCC] <- AnzGeneDirektAnnotiertAllUniq[NichtNA]
AnzGeneDirektAnnotiertManUniqCC[IndexElternCC] <- AnzGeneDirektAnnotiertManUniq[NichtNA]
CCAnnotInd <- match(termId(KeyCC),names(AnzAllAnnotations))
AnzAllAnnotationsCC <- KeyCC*0
AnzAllAnnotationsCC[which(!is.na(CCAnnotInd))] <- AnzAllAnnotations[CCAnnotInd[which(!is.na(CCAnnotInd))]]
Pti <-  AnzAllAnnotationsCC/AnzAllAnnotations[which(names(AnzAllAnnotations)=='GO:0005575')] # Annot. an Term / Annot. an Wurzel
CCAnnotManuInd <- match(termId(KeyCC), names(AnzManuAnnotations))
AnzAllManuAnnotationsCC <- KeyCC*0
AnzAllManuAnnotationsCC[which(!is.na(CCAnnotManuInd))] <- AnzManuAnnotations[CCAnnotManuInd[which(!is.na(CCAnnotManuInd))]]
PtiManu <- AnzAllManuAnnotationsCC/AnzManuAnnotations[which(names(AnzManuAnnotations)=='GO:0005575')]
#Minimum ist 1/17839 = 5.605695e-05
#Dementsprechend Maximum InformationContent = -log2(1/17839) = 14.12275. (gerundet 14.12)
#Idee: Evtl sollte man das zur Vergleichbarkeit mit anderen Ontologien in Prozent von Maximum angeben... 
NullIndex <- which(Pti==0)
InformationContentCC <- round(x = -log2(Pti), digits = 2)
InformationContentCC[NullIndex] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = -1.
NullIndexManu <- which(PtiManu==0)
InformationContentCCManu <- round(x = -log2(PtiManu), digits = 2)
InformationContentCCManu[NullIndexManu] <- NaN # Fuer die Terme ohne Annotierung setzen wir IC = NaN.
AlleAnnotierungenCC <- AnzAllAnnotationsCC
TiefeCC <- MinPathL2Root('CC')
LevelCC <- MaxPathL2Root('CC')
HoeheCC <- MinPathL2Leaf('CC')
HoehenLevelCC <- MaxPathL2Leaf('CC')

# erstelle aus Daten Matrix
Data <- cbind(KeyCC, AnzGeneDirektAnnotiertCC, AlleAnnotierungenCC, AnzGeneDirektAnnotiertManCC, AnzGeneDirektAnnotiertAllUniqCC, AnzGeneDirektAnnotiertManUniqCC, AnzAllManuAnnotationsCC, AnzahlNachkommenCC, AnzahlKinderCC, AnzahlVorfahrenCC, AnzahlElternCC, InformationContentCC, InformationContentCCManu, TiefeCC, LevelCC, HoeheCC, HoehenLevelCC)
#Schreibe aus.
suppressWarnings(WriteLRN('GOTermInfosCC', Data = Data, Header = c('GOtermNr', 'AnzGeneDirektAnnotiertCC', 'AlleAnnotierungenUniqCC', 'AnzGeneDirektAnnotiertManCC', 'AnzGeneDirektAnnotiertAllUniqCC', 'AnzGeneDirektAnnotiertManUniqCC', 'AnzAllManuAnnotationsCC', 'AnzahlNachkommenCC', 'AnzahlKinderCC', 'AnzahlVorfahrenCC', 'AnzahlElternCC', 'InformationContentCC', 'InformationContentCCManu', 'MinPathL2RootCC', 'MaxPathL2RootCC', 'MinPathL2LeafCC', 'MaxPathL2LeafCC'), 
Key = KeyCC, DataDefined = c(), OutDirectory = OutDirectory, CommentOrDigits = Comments))
print('16/18 done.')
 
 
 
 
 
##########################################################################################
## AllAnnNCBIsPlusGeneName.names #########################################################
##########################################################################################
# Schreibe alle NCBIs, die zu mindestens einem GO-Term annotiert sind aus. (Plus GeneNames)
# Alle NCBIs finden, die zu mind. einem GO-Term annotiert sind:
allNCBIsUnsort <- unique(unlist(unname(GO2ALLEGS))) # Liste auseinandernehmen und nur die NCBIs behalten als Strings
allNCBIs <- sort(as.numeric(allNCBIsUnsort))
# Descriptions der Gene finden.
y <- org.Hs.eg.db::org.Hs.egGENENAME
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- AnnotationDbi::mappedkeys(y)
# Convert to a list
yy <- as.list(y[mapped_genes])
# GeneNames (Also die Symbols) fuer die NCBIs finden
z <- org.Hs.eg.db::org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- AnnotationDbi::mappedkeys(z)
# Convert to a list
zz <- as.list(z[mapped_genes])
# Ausschreiben:
GeneNames <- unlist(yy[as.character(allNCBIs)])
GeneSymbols <- unlist(zz[as.character(allNCBIs)])
FileName <- 'AllAnnNCBIsPlusGeneName.names'
Names <- GeneSymbols
Key <- allNCBIs
FurtherTexts <- GeneNames
DescriptionHeader <- c('NCBI','GeneNames\tGeneDescription')
Comments <- paste0('All NCBIs (as Key) that are annotated to at least one GO term with corresponding gene names and description.\n# Last update: ', Sys.Date())
suppressWarnings(WriteNAMES(FileName, Names, Key, FurtherTexts, OutDirectory, DescriptionHeader, Comments))
print('17/18 done.')


##########################################################################################
## EvidenceCodes.names #########################################################
##########################################################################################
# Erstelle Datei mit allen Evidence Codes mit denen NCBIs zu GO-Termen annotiert sein koennen.
# Erzeugen dabei numerischen Key, der eindeutig den Evidence-Code-Strings zugeordnet ist.

OccuringEvidenceCodes <- unique(unlist(unname(lapply(GO2ALLEGS, names))))

# Nur IEA ist automatisch curatiert, alle anderen Evidence Codes sind manuell curatiert. Daher
# schreiben wir IEA an erste Stelle und alle anderen Codes dahinter in alphabetischer Reihenfolge.
 
ieaInd <- which(OccuringEvidenceCodes=='IEA')
OccEC <- OccuringEvidenceCodes[-ieaInd]
Names <- c('IEA',sort(OccEC))
Key <- 2^(0:(length(Names)-1))
FileName <- 'EvidenceCodes.names'
OutDirectory <- OutDirectory
DescriptionHeader <- c('Key', 'EvidenceCode')
Comments <- paste0('List of EvidenceCodes that occur in org.Hs.eg.db package with unique Key. \n# Date: ', Sys.Date(), '. \n# From the documentation of org.Hs.eg.db package: \n# The Evidence element contains a code indicating what kind of evidence supports the association of \n# the GO identifier to the Entrez Gene id. Some of the evidence codes in use include: \n# IMP: inferred from mutant phenotype \n# IGI: inferred from genetic interaction \n# IPI: inferred from physical interaction \n# ISS: inferred from sequence similarity \n# IDA: inferred from direct assay \n# IEP: inferred from expression pattern \n# IEA: inferred from electronic annotation \n# TAS: traceable author statement \n# NAS: non-traceable author statement \n# ND: no biological data available \n# IC: inferred by curator \n#' )

# WriteNAMES(FileName, Names,Key,FurtherTexts,OutDirectory, DescriptionHeader, Comments)
suppressWarnings(WriteNAMES(FileName, Names, Key, OutDirectory=OutDirectory, DescriptionHeader=DescriptionHeader, Comments=Comments))
print('18/18 done.')

}# end function updateORAdatabase