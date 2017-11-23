ORAfilename <- function(OutFile, NrOfValidInputGenes, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur,# ABCAnalysis, ImpThreshold
WithRefSet = FALSE){
# Function to complete the OutFile name passed to dbtORA with ORA parameters.

# OutFilePlusParams <- ORAfilename(OutFile, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur)
# # # OutFilePlusParams <- ORAfilename(OutFile, Correction, PvalueThreshold, MinNrOfGenes, OnlyManuCur, ABCAnalysis, ImpThreshold)

# INPUT:
# Most of the parameters given to dbtORA (default values set in dbtORA!):
# OutFile								String;
#												Filename of the output file(s). Will be complemented by the parameters of the ORA.
# NrOfValidInputGenes		String;
#												Number of valid input genes = #(input genes) - #(duplicated and non-annotated genes) 
# Correction						String;
#												Type of correction for mulitple testing of the p-values.
#												'BON' for Bonferroni, 
#												'FDR' for False Discovery Rate,
#												'RAW' if no correction should be done.
# PvalueThreshold				Numeric;
#												P-value threshold. GO-Terms with p-values greater than PvalueThreshold will be ignored.
# MinNrOfGenes					Numeric;
#												Minimum number of genes annotated to one Term that is accepted. GO-Terms with more than
#												MinNrOfGenes genes will be considered in calculation.
# OnlyManuCur						Boolean;
#												Set TRUE if only manually curated gene annotations should be considered.
# # # ABCAnalysis							Boolean;
# # #													Do an ABC analysis to curt the vector of GO terms used to draw the DAG to group A terms.
# # # ImpThreshold						Numeric;
# # #													Importance Threshold. If the calculated importance of a term is smaller than ImpThreshold,
# # #													the term won't occur in the DAG.
# WithRefSet						Boolean;
#												Set TRUE if RefSet given. Default: FALSE.			

# OUTPUT:
# OutFilePlusParams			String;
#												The complemented OutFile.

# AUTHOR:
# CL, 03.03.2016

# USES: 
# [1] "function"  "fileparts" "gsub"      "paste"     "if"        "is.null"  
# [7] "return"   


# Wir muessen die Dateiendung entfernen, sofern gegeben, Parameter anfuegen und Dateiendung wieder anhaengen.
OutFileParts <- fileparts(OutFile)
OutFileWithoutExt <- gsub(OutFileParts$ext, '', OutFile)
OutFileWithoutExt <- paste0(OutFileWithoutExt, 'Genes', NrOfValidInputGenes)
OutFilePlusParams <- paste(OutFileWithoutExt, Correction, PvalueThreshold, MinNrOfGenes, sep='_')
if(OnlyManuCur){
	OutFilePlusParams <- paste(OutFilePlusParams, 'MANU', sep = '_')
}else{
	OutFilePlusParams <- paste(OutFilePlusParams, 'ALL', sep = '_')
}# end if(OnlyManuCur)
# if(ABCAnalysis){
# 	OutFilePlusParams <- paste(OutFilePlusParams, 'ABC', sep = '_')
# }
# if(ImpThreshold != 0){
# 	OutFilePlusParams <- paste(OutFilePlusParams, 'ImpThresh', sep = '_'
# }
if(WithRefSet){
	OutFilePlusParams <- paste(OutFilePlusParams, 'RefSet', sep= '_')
}

return(OutFilePlusParams)
}# end function ORAfilename
