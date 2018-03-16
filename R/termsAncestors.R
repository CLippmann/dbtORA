termsAncestors <- function(GOTermNr, OntologyNr){
# Function returns vector of all ancestors in ontology for given GO term numbers.
#
# Ancestors  <- termsAncestors(GOTermNr, OntologyNr) 
#
# NOTA was ist mit den Relationstypen passiert z.B regulates? 
# eigentlich folgen wir nur part-of und is-a
#
# INPUT:
# GOTermNr			Numeric; 
#								Vector of GO term numbers.
# OntologyNr		Numeric;
#								To select the ontology. One of: 
#								1 for biological process, 
#								2 for molecular function or 
#								4 for cellular component
#
# OUTPUT:
# List of 3:
#	Ancestors								Numeric; Unique and by number sorted GO terms that are ancestors of input GO terms.
# TermsWithoutAncestors		Numeric; Vector of those GO terms for that no ancestors could be found.				
#	GOterms2GOtermsMatrix 	Numeric; Sparse matrix of GO terms and ancestors.

# AUTHOR:
# CL 21.03.2016 (Vorlage von ALU, CL & MT: TermAncestors)

# USES: 
# functions:
 # [1] "function"         "if"               "length"           "stop"            
 # [5] "all"              "is.numeric"       "switch"           "as.character"    
 # [9] "ReadSparseMatrix" "system.file"      "termNr"           "match"           
# [13] "which"            "is.na"            "sort"             "unique"          
# [17] "return"           "list"  


# Dateneingabe ueberpruefen:
if(length(GOTermNr) <= 0){stop('termsAncestors: A non-empty vector of GO term numbers is required as input! Function stops.')}
if(!all(is.numeric(GOTermNr))){stop('termsAncestors: GO term numbers are not numeric! Function stops.')}
if(!is.numeric(OntologyNr)){stop('termsAncestors: OntologyNr has to be numeric: 1 (=BP), 2 (=MF) or 4 (=CC)! Function stops.')}

	
# Infos ueber Terme aus GOdata einlesen:
switch(as.character(OntologyNr),
	'1' =	{	# biological process = BP
					AncSparseMatrix <- ReadSparseMatrix('BPAncestorsSparseMatrix', system.file('extdata',package='ORA'))
					#GOroot <- 8150
				},
	'2' = { # molecular function = MF
					AncSparseMatrix <- ReadSparseMatrix('MFAncestorsSparseMatrix', system.file('extdata',package='ORA'))
					#GOroot <- 3674
				},
	'4' = { # cellular component = CC
					AncSparseMatrix <- ReadSparseMatrix('CCAncestorsSparseMatrix', system.file('extdata',package='ORA'))
					#GOroot <- 5575
				},
	stop('termsAncestors: OntologyNr has to be 1 (=BP), 2 (=MF) or 4 (=CC)! Function stops.')
)# end switch(OntologyNr)
  
# vorlaeufige Liste aller Vorfahren ungeordnet und mit Duplikaten erstellen:
AncSparseNames <- termNr(AncSparseMatrix$DimNames[[1]]) # Dimnames sind GOtermIDs - umwandeln in Nummern
GOInd <- match(GOTermNr, AncSparseNames) # Index nach dem wir ZeilenInd durchsuchen muessen -> umrechnen der GOtermNr in die Indizes, die auch in ZeilenInd und SpaltenInd verwendet werden
NaInd <- which(is.na(GOInd)) # Terme, die nicht in der Ontologie sind oder nicht existieren (keine Annotationen haben und somit fuer uns uniteressant sind.)
NotNaInd <- which(!is.na(GOInd)) # Terme, die existieren
TermeOhneEltern <- GOTermNr[NaInd] # da wo bei match NA rauskommt, koennen wir keine Eltern finden
# GOIndOhneNa <- GOInd[-NaInd] # aber dann stehn die Terme nicht mehr an der richtigen stelle... iwie merken... 

# Welche ZeilenInds stimmen mit GOInds (ohne NAs) ueberein?
ZeilenIndInGOIndMitNAs <- match(AncSparseMatrix$ZeilenInd, GOInd[NotNaInd]) # kann keine NAs geben.
NotNaIndZeilenIndInGOInd <- which(!is.na(ZeilenIndInGOIndMitNAs))
# ohne NAs:
ZeilenIndInGOInd <- ZeilenIndInGOIndMitNAs[NotNaIndZeilenIndInGOInd]
# zugehoerigen SpaltenInd:
SpaltenIndInGOInd <- AncSparseMatrix$SpaltenInd[NotNaIndZeilenIndInGOInd] # hier sind dann die ganzen Ancestors drin
# InputGOtermNrOhneNA <- unique(AncSparseNames[GOInd[NotNaInd]])
Ancestors4InputGOtermNr <- sort(unique(AncSparseNames[SpaltenIndInGOInd]))# GOtermNrs der Ancestors

# if(length(TermeOhneEltern)>0){ 
	# print(paste0('termsAncestors(): No ancestors found for the following GO-terms (represented by GO-term numbers): ', TermeOhneEltern,'.'))
# }# end if(length(TermeOhneEltern)>0)

return(list(Ancestors = Ancestors4InputGOtermNr, TermsWithoutAncestors = TermeOhneEltern))
}# end function termsAncestors

