ask2loadFile <-function(extension, InDirectory = ''){
# V <-ask2loadFile(".lrn")
# FileName = V$FileName
# InDirectory = V$InDirectory
# Benutzerdialog zum Einstellen eines Directory und zum Vorbereiten des
# Lesens des entsprechenden DataIO Dateityps
#
# INPUT
# extension             Angabe welcher dateityp gelesen werden soll Beispiel:  extension ='lrn'
# InDirectory           Startordner in welchem gesucht werden soll
#
# OUTPUT
# FileName          Dateiname der vom Benutzer ausgeaaehlt wurde
# InDirectory         Directory welches vom Benutzer bestimmt wurde
# author MT 12/2015
# Nota: Im Gegensatz zu Matlab wird kein DialogFilename benoetigt

# Edit: CL 24.01.18

if(length(extension)==0){stop('ask2loadFile: No valid extension!')}

extension = gsub('\\.','', extension) #erstmal alle Punkte entfernen
extension = gsub('\\*','', extension) #alle Sternchen entfernen
extension = paste0('*.',extension)		#ganz ordentlich nur ein Sternchen und einen Punkt vorne anfuegen

if(length(extension)==1){
	if(length(InDirectory)>0){InDirectory = paste0(dirname(InDirectory), extension)}
	res=NULL
	tryCatch({
		if(Sys.info()["sysname"]=="Windows"){
			filter=c(extension,extension)
			pathLRN <- utils::choose.files(default=InDirectory, caption=paste("Choose",extension,"File"), multi=FALSE, filters=filter)
		}else{
			print(paste("Please choose a ",extension,"-File."))
			pathLRN <- file.choose()
		}# end if Systemabhaengige File-Abfrage

		if(length(pathLRN) == 0) return(NULL) #Wenn der User seine Eingabe abbricht NULL zurueckgeben
		
		if(!(paste0("*.", tools::file_ext(pathLRN)) %in% extension)){
			warning(paste("You did not select a ",extension,"-File!"))
			return(NULL)
		}
		
		res = list(FileName = basename(pathLRN), InDirectory = dirname(pathLRN))
	}, error = function(ex){
			print('Error while loading data.')
			res=NULL
		}
	)# end tryCatch
	return(res)
	
}else{ # wir haben also mehr als eine Extension uebergeben bekommen
	res=NULL
	extensions <- paste(extension, collapse = ';')
	if(length(InDirectory)>0){InDirectory = paste0(dirname(InDirectory), "/*.*")}
	tryCatch({
		if(Sys.info()["sysname"]=="Windows"){
			filter = c(extensions, extensions)
			pathLRN <- utils::choose.files(default = InDirectory, caption=paste("Choose",extensions,"File"), multi=FALSE, filters=filter)
		}else{
			print(paste("Please choose a ",extensions,"-File."))
			pathLRN <- file.choose()
		}# end if Systemabhaengige File-Abfrage

		if(length(pathLRN) == 0) return(NULL) #Wenn der User seine Eingabe abbricht NULL zurueckgeben

		if(!(paste0("*.",tools::file_ext(pathLRN)) %in% extension)){
			warning(paste("You did not select a",extension,"File!"))
			return(NULL)
		}
		
		res=list(FileName = basename(pathLRN), InDirectory = dirname(pathLRN))
	}, error = function(ex){
			print('Error while loading data.')
			res=NULL
		}
	)# end tryCatch
	return(res)

}# end if mehr als eine Extension gegeben

}# end function ask2loadFile
