checkORAparameters <- function(InFileWithExt, InFileDirectory, RefSetFileWithExt, RefSetDirectory, OutFile, OutFileDirectory, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, drawDAG, MarkDetails, MarkHeadlines, # ABCAnalysis, ImpThreshold,
PlotExt){
# Function to check if the parameters passed to dbtORA are correct.

# checkORAparameters(InFileWithExt, InFileDirectory, RefSetFileWithExt, RefSetDirectory, OutFile, OutFileDirectory, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, drawDAG, MarkDetails, MarkHeadlines, PlotExt)
# # # checkORAparameters(InFileWithExt, InFileDirectory, RefSetFileWithExt, RefSetDirectory, OutFile, OutFileDirectory, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, drawDAG, MarkDetails, MarkHeadlines, ABCAnalysis, ImpThreshold, PlotExt)

# INPUT:
# All parameters of dbtORA (default values set in dbtORA):
# InFileWithExt				String; 
#											Filename with extension where NCBIs are keys (for *.names and *.lrn files)
#											or the only column (for *.txt files).
# InFileDirectory			String;
#											Directory where InFileWithExt can be found.
# RefSetFileWithExt		String;
#											Filename with extension where NCBIs are keys (for *.names and *.lrn files)
#											or the only column (for *.txt files). NCBIs will be used as reference set.
# RefSetDirectory			String;
#											Directory where RefSetFileWithExt with reference NCBIs can be found.
# OutFile							String;
#											Filename of the output file(s). Will be complemented by the parameters of the ORA.
# OutFileDirectory		String;
#											Directory where results of ORA and DAGs will be saved.
# Correction					String;
#											Type of correction for mulitple testing of the p-values.
#											'BON' for Bonferroni, 
#											'FDR' for False Discovery Rate,
#											'RAW' if no correction should be done.
# PvalueThreshold			Numeric;
#											P-value threshold. GO-Terms with p-values greater than PvalueThreshold will be ignored.
# MinNrOfGenes				Numeric;
#											Minimum number of genes annotated to one Term that is accepted. GO-Terms with more than
#											MinNrOfGenes annotated genes will be considered in calculation.
# OnlyManuCur					Boolean;
#											Set TRUE if only manually curated gene annotations should be considered.
# drawDAG							Boolean;
#											Set TRUE if directed acyclic graphs (DAGs) should be drawn.
#											If drawDAG is set to FALSE, the parameters MarkDetails, MarkHeadlines and CurtDAG will be
#											ignored.
# MarkDetails					Boolean;
#											Set TRUE if details of the DAG should be marked in blue colour.
# MarkHeadlines				Boolean;
#											Set TRUE if headlines of the DAG should be marked in yellow colour.
# # # # ABCAnalysis					Boolean;
# # # #											Do an ABC analysis to curt the vector of GO terms used to draw the DAG to group A terms.
# # # # ImpThreshold				Numeric;
# # # #											Importance Threshold. If the calculated importance of a term is smaller than ImpThreshold,
# # # #											the term won't occur in the DAG.
# PlotExt							String;
#											Extension of the plotfile showing the DAG. One of 'pdf' , 'eps' or 'png'.

# AUTHOR:
# CL, 03.03.2016

# USES:
# functions:
 # [1] "function"           "c"                  "is.character"       "is.null"            "is.numeric"        
 # [6] "is.logical"         "rep"                "sum"                "length"             "if"                
# [11] "which"              "stop"               "paste"              "names"              "file_ext"          
# [16] "basename"           "file_path_sans_ext" "dirname"            "dir.exists"         "file.exists"       
# [21] "paste0"             

##################################################
# Datentypen der einzelnen Parameter ueberpruefen:
Datentyp <- c(
		InFileWithExt = is.character(InFileWithExt),
		InFileDirectory = is.character(InFileDirectory),
		RefSetFileWithExt = is.character(RefSetFileWithExt)|is.null(RefSetFileWithExt),
		RefSetDirectory = is.character(RefSetDirectory),
		OutFile = is.character(OutFile),
		OutFileDirectory = is.character(OutFileDirectory),
		Correction = is.character(Correction),
		PvalueThreshold = is.numeric(PvalueThreshold),
		MinNrOfGenes = is.numeric(MinNrOfGenes),
		OnlyManuCur = is.logical(OnlyManuCur),
		drawDAG = is.logical(drawDAG),
		MarkDetails = is.logical(MarkDetails),
		MarkHeadlines = is.logical(MarkHeadlines),
		# ABCAnalysis = is.logical(ABCAnalysis),
		# ImpThreshold = is.numeric(ImpThreshold),
		PlotExt = is.character(PlotExt)
	)#end Datentyp
KorrekterTyp <- c(rep('character', 7), rep('numeric', 2), rep('logical', 4), 'character')
DatentypOK <- sum(Datentyp)==length(Datentyp) # Alle TRUE?
if(!DatentypOK){ # Wenn nicht alle TRUE, gibt die Parameter, die den falschen Typ haben aus.
	FALSEparam <- which(Datentyp == FALSE)
	stop('checkORAparameters: Parameter(s) ',paste(names(Datentyp)[FALSEparam], collapse = ', '),' do not have correct type. Has to be of type "', KorrekterTyp[FALSEparam],'". Function stops.')
}

######################################################################################################
# Jetzt die einzelnen Parameter durchgehen und spezifische Eigenschaften abfragen:									 #
######################################################################################################

######################################################################################################
# Ist InFileWithExt ein filename mit einer der zugelassenen Extensions im Verzeichnis InFileDirectory? 
InFileWithExtExt <- file_ext(InFileWithExt)
InFileWithExtName <- basename(tools::file_path_sans_ext(InFileWithExt))
InFileWithExtDir <- dirname(InFileWithExt)
# Extension von InFileWithExt ok?
if(!(InFileWithExtExt %in% c('names', 'lrn', 'txt'))){
	stop('checkORAparameters: Extension of InFileWithExt is not correct. Has to be one of "lrn", "names" or "txt". Function stops.')
}
# Verzeichnispfad InFileDirectory existent und bei InFileWithExt kein Pfad angegeben?
if(!dir.exists(InFileDirectory)){
	stop('checkORAparameters: InFileDirectory does not exist. Function stops.')
}
if(InFileWithExtDir != '.'){
	stop('checkORAparameters: InFileWithExt may not contain complete path. Function stops.')
}
# File InFileWithExt im Verzeichnis InFileDirectory existent?
FileOK <- file.exists(paste0(InFileDirectory, '/', InFileWithExtName,'.', InFileWithExtExt))
if(!FileOK){
	stop('checkORAparameters: There is no file called ', InFileWithExtName,'.', InFileWithExtExt,' in InFileDirectory. Function stops.')
}

######################################################################################################
# Ist RefSetFileWithExt (falls gegeben) ein filename mit einer der zugelassenen Extensions 
# im Verzeichnis RefSetDirectory?
if(!is.null(RefSetFileWithExt)){
	RefSetFileWithExtExt <- file_ext(InFileWithExt)
	RefSetFileWithExtName <- basename(tools::file_path_sans_ext(InFileWithExt))
	RefSetFileWithExtDir <- dirname(InFileWithExt)
		# Extension von RefSetFileWithExt ok?
	if(!(RefSetFileWithExtExt %in% c('names', 'lrn', 'txt'))){
		stop('checkORAparameters: Extension of RefSetFileWithExt is not correct. Has to be one of "lrn", "names" or "txt". Function stops.')
	}
	# Verzeichnispfad RefSetDirectory existent und bei RefSetFileWithExt kein Pfad angegeben?
	if(!dir.exists(RefSetDirectory)){
		stop('checkORAparameters: RefSetDirectory does not exist. Function stops.')
	}
	if(RefSetFileWithExtDir != '.'){
		stop('checkORAparameters: RefSetFileWithExt may not contain complete path. Function stops.')
	}
	# File RefSetFileWithExt im Verzeichnis RefSetDirectory existent?
	FileOK <- file.exists(paste0(RefSetDirectory, '/', RefSetFileWithExtName, '.', RefSetFileWithExtExt))
	if(!FileOK){
		stop('checkORAparameters: There is no file called ', RefSetFileWithExtName, '.', RefSetFileWithExtExt,' in RefSetDirectory. Function stops.')
	}
}# end if(!is.null(RefSetFileWithExt))

######################################################################################################
# Existiert OutFileDirectory?
if(!dir.exists(OutFileDirectory)){
	stop('checkORAparameters: OutFileDirectory does not exist. Function stops.')
}

######################################################################################################
# Ist Correction "BON", "FDR" oder "RAW" oder etwas unerlaubtes?
if(!(Correction %in% c('BON', 'FDR', 'RAW'))){
	stop('checkORAparameters: Correction has to be one of "BON", "FDR" or "RAW". Function stops.')
}

######################################################################################################
# Ist P-Value Threshold im Intervall von 0 bis 1?
if(!(0 <= PvalueThreshold & PvalueThreshold <= 1)){
	stop('checkORAparameters: P-value threshold has to be in the interval [0;1]. Function stops.')
}

######################################################################################################
# Ist MinNrOfGenes groesser gleich Null?
if(MinNrOfGenes < 0){
	stop('checkORAparameters: MinNrOfGenes has to be a positive number. Function stops.')
}

######################################################################################################
# Ist Importance Threshold zwischen 0 und 100?
# if(!(0 <= ImpThreshold & ImpThreshold <= 100)){
	# stop('checkORAparameters: Importance threshold has to be in the interval [0;100]. Function stops.')
# }

######################################################################################################
# Ist die Extension, die fuer den Plot benutzt werden soll, zulaessig?
if(!(PlotExt %in% c('png', 'eps', 'pdf'))){
	stop('checkORAparameters: PlotExt has to be one of "png", "eps" or "pdf". Function stops.')
}

}#end function checkORAparameters