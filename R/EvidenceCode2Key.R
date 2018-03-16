EvidenceCode2Key <- function(EvidenceCode){
# Function to get the numeric Key for character EvidenceCodes.
# EvidenceKey <- EvidenceCode2Key(EvidenceCode)
#
# INPUT:
# EvidenceCode	String vector[1:n] containing EvidenceCodes. 
#								Known Codes:	"ND"  "IDA" "IEA" "TAS" "NAS" "IMP" "IPI" 
#															"IBA" "ISS" "IEP" "IC"  "IGI" "EXP"
# OUTPUT:
# EvidenceKey		Numeric vector[1:n] containing corresponding Keys according to
#								EvidenceCodes.names.
#
# USES:
# ReadNAMES
#
# AUTHOR:
# CL 28.10.2015

Evi <- ReadNAMES('EvidenceCodes.names', system.file('extdata',package='ORA'))
MatchInd <- match(EvidenceCode, Evi$Names)
EvidenceKey <- Evi$Key[MatchInd]
if(any(is.na(EvidenceKey))){
	print('evidence2Key: One or more input Evidence Codes are unknown. Output contains NAs.')
}
return(EvidenceKey)
}#end function