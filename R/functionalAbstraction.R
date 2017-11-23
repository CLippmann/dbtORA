functionalAbstraction <- function(ORAlrnFilename, ORAdirectory){
# Function to do a naive functional abstraction, i.e. taking the headlines 
# of ORA results as functional areas or if there are too many, subsuming them
# by taking their parents as functional areas.

# functionalAbstraction(ORAlrnFile, ORAdirectory)

# INPUT:
# OPTIONAL:
# ORAlrnFilename	String; Filename of the *.lrn-file from ORA output.
# ORAdirectory		String; Directory where ORAlrnFilename can be found.

# OUTPUT:
# Plots of sub-DAGs with the functional areas as roots and a *.names-file 
# containing the functional areas.

if(missing(ORAlrnFilename)){ 
	# interaktiv abfragen, wo die Datei ist.
	# Frage Werte fuer ORAlrnFilename und ORAdirectory ab.
	print('Please select *lrn-file containing the results of ORA.')
	flush.console() # print wird sofort ausgegeben.
	File <- ask2loadFile("lrn")
	ORAlrnFilename <- File$FileName
	ORAdirectory <- File$InDirectory
	print(paste0(ORAlrnFilename, ' selected.'))
}
ORAlrn <- ReadLRN(ORAlrnFilename, ORAdirectory)
# > str(ORAlrn)
# List of 5
 # $ Key        : num [1:140] 6796 6793 30182 50877 42325 ...
 # $ Data       : num [1:140, 1:16] 1 1 1 1 1 1 1 1 1 1 ...
 # $ Header     : chr [1:16] "OntologyNr" "NrOfGenesInUniverse" "NrOfGenesInSample" "NrOfAnnotationsInTerm" ...
 # $ DataDefined: num [1:17] 9 1 1 1 1 1 1 1 1 1 ...
 # $ Comments   : chr " WriteORAresults: ORAresults$LRNresults.    Parameters: ORA,FDR,0.05,2,ALL    Original input file: 01HSANgene11.names"
ORAfile <- fileparts(ORAlrnFilename)$name # Nur der Name ohne Extension
# Jetzt noch die Zahlen am Ende loswerden:
# tmp <- attr(gregexpr('\\(?[0-9,.]+', ORAfile)[[1]], "match.length") # liefert die Laengen der vorkommenden Zahlen 
# # brauchen nur den letzten Eintrag in tmp
# ORAfile <- strtrim(ORAfile, nchar(ORAfile)-tmp[length(tmp)])
ORAfile2 <- strtrim(ORAfile, nchar(ORAfile)-nchar(dim(ORAlrn$Data)[1]))
ORAnames <- ReadNAMES(ORAfile, ORAdirectory)
# > str(ORAnames)
# List of 4
 # $ Key         : num [1:140] 6796 6793 30182 50877 42325 ...
 # $ Names       : chr [1:140] "phosphate-containing compound metabolic process" "phosphorus metabolic process" "neuron differentiation" "neurological system process" ...
 # $ FurtherTexts: chr [1:140, 1] "GO:0006796" "GO:0006793" "GO:0030182" "GO:0050877" ...
 # $ Comments    : chr " WriteORAresults: ORAresults$NAMESresults.    Parameters: ORA,FDR,0.05,2,ALL    Original input file: 01HSANgene11.names GOtermN"| __truncated__
if(!TheSameKey(ORAlrn$Key, ORAnames$Key)){stop('functionalAbstraction: ORA *.lrn-file and ORA *.names-file do not have the same key. Function stops.')}
# Headlines rausfinden:
HeadlineInd <- which(ORAlrn$Header == "IsHeadline")
# FAdata <- ORAlrn$Key[as.logical(ORAlrn$Data[,HeadlineInd]),]
FAKey <- ORAnames$Key[as.logical(ORAlrn$Data[,HeadlineInd])]
FANames <- ORAnames$Names[as.logical(ORAlrn$Data[,HeadlineInd])]
FAFurtherTexts <- ORAnames$FurtherTexts[as.logical(ORAlrn$Data[,HeadlineInd])]
FAFileName <- paste0(ORAfile2,'_NaiveFA', length(FAKey))
# WriteNAMES(FileName = ,Names = ,Key = ,FurtherTexts = ,OutDirectory = ,DescriptionHeader = ,Comments = )
WriteNAMES(FAFileName, FANames, FAKey, FAFurtherTexts, ORAdirectory, DescriptionHeader = c('GOtermNr', 'GOtermDescription', 'GOtermID'), Comments = 'Naive functional areas = headlines from ORA.')

AdjFileName <- list.files(ORAdirectory, pattern = paste0(ORAfile2, 'GOterms2GOterms')) # Mit den gegebenen Parametern aus ORA kann es nur eine Adj geben - wir kennen nur die Dimension nicht und muessen deshalb so suchen!
# drawFunctionalAreas(FAfileName = , FAdir = , ORAfileName = , AdjFileName = , ORAdir = , PlotDir = , Ontology = , DrawRoot2FAs = , MarkHeadlines = , MarkDetails = , PlotFileExt = )
drawFunctionalAreas(FAFileName, ORAdirectory, ORAlrnFilename, AdjFileName)

### ganz naiv. Viel zu viele FAs... muessten wenigstens noch iwie headlines zusammenfassen, die nur
# headlines als eltern haben. oder so.
}