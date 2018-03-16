dbtORA <- function(InFileWithExt, PvalueThreshold = 0.05, Correction = 'BON', OnlyManuCur = TRUE, MinNrOfGenes = 2, InFileDirectory = getwd(), OutFile = InFileWithExt, OutFileDirectory = InFileDirectory, RefSetFileWithExt = NULL, RefSetDirectory = InFileDirectory, drawDAG = TRUE, MarkDetails = TRUE, MarkHeadlines = TRUE,  #ABCAnalysis = FALSE, ImpThreshold = 0, 
PlotExt = 'png'){ #, SuggestParams = FALSE){
# Function to do an Overrepresentation Analysis including the drawing of DAGs.

# V <- dbtORA(InFileWithExt, PvalueThreshold = 0.05, Correction = 'BON', OnlyManuCur = TRUE, MinNrOfGenes = 2, InFileDirectory = getwd(), OutFile = InFileWithExt, OutFileDirectory = InFileDirectory, RefSetFileWithExt = NULL, RefSetDirectory = InFileDirectory, drawDAG = TRUE, MarkDetails = TRUE, MarkHeadlines = TRUE, PlotExt = 'png')
# # # V <- dbtORA(InFileWithExt, PvalueThreshold = 0.05, Correction = 'BON', OnlyManuCur = TRUE, MinNrOfGenes = 2, InFileDirectory = getwd(), OutFile = InFileWithExt, OutFileDirectory = InFileDirectory, RefSetFileWithExt = NULL, RefSetDirectory = InFileDirectory, drawDAG = TRUE, MarkDetails = TRUE, MarkHeadlines = TRUE, ABCAnalysis = FALSE, ImpThreshold = 0, PlotExt = 'png')

# INPUT:
# InFileWithExt			String; 
#										Filename with extension where NCBIs are keys (for *.names and *.lrn files)
#										or the only column (for *.txt files).

# OPTIONAL:
# PvalueThreshold		Numeric; Default: 0.05
#										P-value threshold. GO terms with p-values greater than PvalueThreshold will be ignored.
# Correction				String; Default: 'BON'
#										Type of correction for mulitple testing of the p-values.
#										'BON' for Bonferroni, 
#										'FDR' for False Discovery Rate,
#										'RAW' if no correction should be done.
# OnlyManuCur				Boolean; Default: TRUE
#										Set TRUE if only manually curated gene annotations should be considered.
# MinNrOfGenes			Numeric; Default: 2
#										Minimum number of genes annotated to one Term that is accepted. Only GO terms with more than
#										MinNrOfGenes annotated genes will be considered in calculation.
# InFileDirectory		String; Default: current directory getwd().
#										Directory where InFileWithExt can be found.
# OutFile						String; Default: InFileWithExt (extension will be adjusted)
#										Filename of the output file(s). Will be complemented by the parameters of the ORA.
# OutFileDirectory	String; Default: InFileDirectory
#										Directory where results of ORA and DAGs will be saved.
# RefSetFileWithExt	String; Default: NULL
#										Filename with extension where NCBIs are keys (for *.names and *.lrn files)
#										or the only column (for *.txt files). NCBIs will be used as reference set.
# RefSetDirectory		String; Default: InFileDirectory
#										Directory where RefSetFileWithExt with reference NCBIs can be found.
# drawDAG						Boolean; Default: TRUE
#										Set TRUE if DAGs should be drawn.
#										If drawDAG is set to FALSE, the parameters MarkDetails, MarkHeadlines and CurtDAG will be
#										ignored.
# MarkDetails				Boolean; Default: TRUE
#										Set TRUE if details of the DAG should be marked in blue colour.
# MarkHeadlines			Boolean; Default: TRUE
#										Set TRUE if headlines of the DAG should be marked in yellow colour.
# # # ABCAnalysis				Boolean; Default: FALSE
# # #										Do an ABC analysis to curt the vector of GO terms used to draw the DAG to group A terms.
# # # ImpThreshold			Numeric; Default: 0
# # #										Importance Threshold. If the calculated importance of a term is smaller than ImpThreshold,
# # #										the term won't occur in the DAG.
# PlotExt						String; Default: 'png' 
#										Extension of the plotfile showing the DAG. One of 'pdf', 'eps' or 'png'.
# # # SuggestParams			Boolean; Default: FALSE
# # #										Set TRUE if suggestions for PvalueThreshold and MinNrOfGenes should be calculated and used for ORA.
# # #										Note: If TRUE, ORA will be calculated twice. This will take some time.

# # # OUTPUT:
# # # *.lrn, *.names, AdjMatrix GO2GO, Matrix Genes2Terms und DAGs (falls drawDAG = TRUE) fuer die nicht interaktive Version.

# AUTHOR:
# CL (uses part of program that was partly written by MT)
# February 2016, version 3

# USES:
# packages: matrixStats, ggm, Matrix, GO.db, tools, utils
# functions:
# [1]  "function"           "getwd"              "c"                  "for"                "suppressMessages"  
 # [6] "requireNamespace"   "if"                 "missing"            "print"              "keys"              
# [11] "flush.console"      "ask2loadFile"       "paste0"             "readline"           "checkORAparameters"
# [16] "switch"             "file_ext"           "ReadNAMES"          "ReadLRN"            "ReadTXT"           
# [21] "unname"             "stop"               "all"                "is.numeric"         "sort"              
# [26] "seq_along"          "warning"            "system.file"        "is.null"            "ORA"               
# [31] "setwd"              "dim"                "length"                        
# [36] "ORAfilename"        "WriteORAresults"    "drawORA"  


# TO DO:
# CurtDAG Filter schreiben. 
# ImpthresholdFilter


# # Lade benoetigte packages (unterdruecke Warnmeldungen, Statusberichte und Infos):
# importantPackages = c("GO.db",  "AnnotationDbi",  "tools", "utils", "Matrix")#"Biobase", "IRanges", "S4Vectors","BiocGenerics",
#for(i in importantPackages) suppressMessages(requireNamespace(i, character.only = T, quietly = T, warn.conflicts = F))

if(missing(InFileWithExt)){ 
	# interaktive Version aufrufen.
	# Frage Werte fuer InFileWithExt, InFileDirectory, RefSetFileWithExt und RefSetDirectory ab.
	# Alle anderen Parameter werden mit Defaultwert verwendet, wenn durch User nicht beim Aufruf spezifiziert.
	print('Please select file containing NCBIs as keys (for *.names and *.lrn files) or NCBIs as the only column (for *.txt files).')
	flush.console() # print wird sofort ausgegeben.
	File <- ask2loadFile(c("names", "lrn", "txt"))
	InFileWithExt <- File$FileName
	InFileDirectory <- File$InDirectory
	print(paste0(InFileWithExt, ' selected.'))
	
  n <- readline("Do you want to use a reference set? y/n  ")
  if(n == 'y' || n == 'Y'){
		print('Please select file containing NCBIs as keys (for *.names and *.lrn files) or NCBIs as the only column (for *.txt files).')
		flush.console() # print wird sofort ausgegeben.
		Refset <- ask2loadFile(c("names", "lrn", "txt"))
		RefSetFileWithExt <- Refset$FileName
		RefSetDirectory <- Refset$InDirectory
		print(paste0(RefSetFileWithExt, ' selected.'))
	}else{ # keine Eingabe oder Mist --> kein Refset
		RefSetFileWithExt <- NULL
		print('ORA without reference set.')
	}
}	

# Sollen Vorschlaege fuer P-value threshold und MinNrOfGenes gemacht werden?
# if(missing(SuggestParams)){
	# m <- readline('Would you like to get suggestions for p-value threshold and minimum number of genes that should be annotated to a GOterm to be considered in DAG? y/n   ')
	# if(m == 'y'||'Y'){
		# SuggParams <- TRUE
	# }else{
		# SuggParams <- FALSE
	# }# end if SuggParams = TRUE
# }

# nichtinteraktive Version:
# Ueberpruefen, ob die Parameter, die der ORA uebergeben wurden, korrekt sind: 
# Falls nicht, bricht Programm ab!
checkORAparameters(InFileWithExt, InFileDirectory, RefSetFileWithExt, RefSetDirectory, OutFile, OutFileDirectory, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, drawDAG, MarkDetails, MarkHeadlines, PlotExt)

# Daten (NCBIs) einlesen:
switch(file_ext(InFileWithExt), 
	'names'	= {	NAMES <- ReadNAMES(InFileWithExt, InFileDirectory)
								NCBIs <- NAMES$Key},
	'lrn'		= {	LRN <- ReadLRN(InFileWithExt, InFileDirectory)
								NCBIs <- LRN$Key},
	'txt'		= {	TXT <- ReadTXT(InFileWithExt, InFileDirectory)
								NCBIs <- unname(TXT$Data[,1])},
	stop('Invalid input format for InFileWithExt!')
)#end switch

# Soweit es geht ueberpruefen, ob die NCBIs auch wirklich NCBIs sind:
if(!all(is.numeric(NCBIs))){
	stop('Non numeric input for NCBIs found in InFileWithExt. Key respectively first column must contain the NCBIs. Function stops.')
}
if(all(sort(NCBIs) == seq_along(NCBIs))){
	warning('Make sure your InFileWithExt contains NCBIs as key respectively as first column. Currently NCBIs are a sequence form 1 to number of NCBIs.')
}

# GOAall einlesen, weil wir sonst nicht wissen welche NCBIs bekannt sind:
GOAall <- ReadLRN('GOAall.lrn',system.file('extdata',package='ORA'))
AllNCBIs <- unname(GOAall$Data[,1])

# Jetzt das Reference Set einlesen, falls gegeben und ORA durchfuehren:
if(is.null(RefSetFileWithExt)){ # RefSet nicht gegeben:
	# Die Berechnung der ORA durchfuehren:
	ORAresults <- ORA(NCBIs, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, RefSet = NULL, GOAall)
}else{ # RefSet gegeben!
	# RefSet einlesen:
	switch(file_ext(RefSetFileWithExt), 
		'.names'	= {	NAMES <- ReadNAMES(RefSetFileWithExt, RefSetDirectory)
									RefSet <- NAMES$Key},
		'.lrn'		= {	LRN <- ReadLRN(RefSetFileWithExt, RefSetDirectory)
									RefSet <- LRN$Key},
		'.txt'		= {	TXT <- ReadTXT(RefSetFileWithExt, RefSetDirectory)
									RefSet <- unname(TXT$Data[,1])},
		stop('Invalid input format for RefSetFileWithExt!')
	)#end switch
	# Soweit es geht ueberpruefen, ob die NCBIs im RefSet auch wirklich NCBIs sind:
	if(!all(is.numeric(RefSet))){
		stop("Non numeric input for reference set's NCBIs found in InFileWithExt. Key respectively first column must contain the NCBIs. Function stops.")
	}
	if(all(sort(RefSet) == seq_along(RefSet))){
		warning('Make sure your RefSetFileWithExt contains NCBIs as key respectively as first column. Currently NCBIs are a sequence form 1 to number of NCBIs.')
	}
	# Die Berechnung der ORA mit RefSet durchfuehren:
	setwd(OutFileDirectory) # damit wir in ORA die doppelten/ungefundenen NCBIs ausschreiben koennen
	ORAresults <- ORA(NCBIs, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, RefSet, GOAall)
}# end if RefSet gegeben oder nicht

# # # # Falls Parameter vorgeschlagen werden sollen, hier berechnen und neue ORA damit durchrechnen:
# # # if(SuggestParams == TRUE){
	# # # V <- suggestORAparams(ORAresults)
	# # # SuggestedPvalThreshold <- V$SuggestedPvalThreshold
	# # # SuggestedMinNrOfGenes <- V$SuggestedMinNrOfGenes
	# # # ## SuggestedRemarkableness <- V$SuggestedRemarkableness
	# # # #IndicesOfRejectedGOterms <- V$IndicesOfRejectedGOterms
	
	# # # print(paste0('Suggested PvalueThreshold: ', SuggestedPvalThreshold, '.'))
	# # # print(paste0('Suggested MinNrOfGenes: ', SuggestedMinNrOfGenes, '.'))
	# # # #print(paste0('By using suggested parameters ', length(IndicesOfRejectedGOterms), ' GOterms are removed from analysis.'))
	# # # print('Calculating ORA with suggested p-value threshold and suggested minimum number of genes that should at least be annotated to one GOterm to be considered in analysis.')
	
	# # # # ORA nochmal mit den neuen Parametern rechnen:
	# # # dbtORA(InFileWithExt = InFileWithExt, PvalueThreshold = SuggestedPvalThreshold, Correction = Correction, OnlyManuCur = OnlyManuCur, MinNrOfGenes = SuggestedMinNrOfGenes, InFileDirectory = InFileDirectory, OutFile = OutFile, OutFileDirectory = OutFileDirectory, RefSetFileWithExt = RefSetFileWithExt, RefSetDirectory = RefSetDirectory, drawDAG = drawDAG, MarkDetails = MarkDetails, MarkHeadlines = MarkHeadlines, PlotExt = PlotExt, SuggestParams = FALSE)
	# # # return()
# # # }# end if SuggestParams TRUE

# Ueberpruefen wie viele Gene wir tatsaechlich zur Berechnung benutzt haben. Doppelte und
# ungefundene werden nicht beruecksichtigt.
# Achtung: Auch "OnlyManuCur" beruecksichtigen - dadurch koennen Gene rausfliegen!
AnzValidInputGenes <- dim(ORAresults$Genes2GOtermsMatrix)[1] -1 # Die Dimension der Genes2Terms-Matrix liefert Anzahl verwendeter Gene

if(AnzValidInputGenes != length(NCBIs)){
	warning(paste0('dbtORA: ', length(NCBIs)-AnzValidInputGenes, ' input gene(s) were not used. There might be duplicates in input genes or some input genes are not annotated to any GO term. For analysis used genes can be found in ', basename(file_path_sans_ext(InFileWithExt)), 'Genes', AnzValidInputGenes, '.names. in ', OutFileDirectory, '.'))
	
	# # Nur benutzte Gene ausschreiben:
	# FileName <- paste0(basename(tools::file_path_sans_ext(InFileWithExt)), 'Genes', AnzValidInputGenes, '.names')
	# Key <- dimnames(ORAresults$Genes2GOtermsMatrix)[[1]][-1]
	# GeneTexts <- NCBI2GeneName(Key)
	# Names <- GeneTexts[[2]]
	# FurtherTexts <- GeneTexts[[1]]
	# OutDirectory <- OutFileDirectory
	# DescriptionHeader <- c('NCBI', 'GeneSymbol', 'GeneName')
	# Comments <- 'Genes from ORA input used for analysis. GeneNames from "AllAnnNCBIsPlusGeneName.names" in system.file('extdata',package='ORA').'
	# WriteNAMES(FileName, Names, Key, FurtherTexts, OutDirectory, DescriptionHeader, Comments)
}# end if duplicated or not found genes exist

# Output-Dateinamen um die Parameter ergaenzen:
if(is.null(RefSetFileWithExt)){ # RefSet nicht gegeben
	OutFilePlusParams <- ORAfilename(OutFile, AnzValidInputGenes, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur)
}else{ # RefSet gegeben
	OutFilePlusParams <- ORAfilename(OutFile, AnzValidInputGenes, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, WithRefSet = TRUE)
}

# Ergebnisse der ORA speichern:
WriteORAresults(OutFilePlusParams, ORAresults, OutFileDirectory, InFileWithExt)

# # # # GOterms2Genes Matrix ausscheiben:
# # # OutMatrixData <- rbind(termNr(RelGOterms),Genes2GOTermsMatrix)  # in die erste Zeile die GOtermNummer
# # # OutMatrixKey  <- c(0,AllAnnNCBIs)            # Schluessel = NCBIs; 1. Zeile ==0
# # # OutMatrixHeader <- c('NCBIs', gsub(' ', '_', termDescription(RelGOterms)))
# # # OutMatrixFilename = paste0(Filename,'GoTerms',AnzahlTerme,'Genes2GOTermMatrix.lrn')
# # # MatrixComment = 'Matrix NCBIs -> GO-terme aus ORA mit NCBIs; Parameter: vgl. Dateiname.\n# '
# # # WriteLRN(OutMatrixFilename, OutMatrixData, OutMatrixHeader, OutMatrixKey, c(), OutDirectory=Directory, MatrixComment)

# Ergebnisse zeichnen, falls gewuenscht:
if(drawDAG){		
	# # # # # Hier GO-Terme rausschmeissen, die je nach Filter nicht gezeichnet werden sollen, oder 
	# # # # # nichts machen, wenn nichts gefiltert werden soll:
	# # # # FilteredORAresults <- filterORAGOterms(ORAresults, ABCAnalysis, ImpThreshold)
	# # # # # Nochmal ausschreiben?
	# # # # # Oder evtl. ne Funktion "plotFilteredORAresults" oder so, wo einfach "filterORAGOterms" aufgerufen wird
	# # # # # und danach plotGOgraph.
	
	# Zeichne jetzt mit allen Infos den DAG und speichere Ergebnisse als Bilder:
	PlotFileWithExt <- paste0(OutFilePlusParams,'.',PlotExt)
	PlotDirectory <- OutFileDirectory
	MarkDetails <- MarkDetails
	MarkHeadlines <- MarkHeadlines
	Overwrite <- TRUE
	drawORA(ORAresults, PlotFileWithExt, PlotDirectory, MarkDetails, MarkHeadlines, Overwrite)

}# end "zeichnen des DAGs" = end if(drawDAG)

}# end function dbtORA
