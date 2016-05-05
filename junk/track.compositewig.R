track.compositewig <-
function( trackfiles , parentname=NULL , hubloc=NULL , range=c(0,3) , plotcolors=rainbow(length(trackfiles)) , tracknames=removeext(trackfiles) ){
	
	if(is.null(parentname)){stop("must define parentname")}
	if(is.null(hubloc)){stop("must define hubloc")}
	trackfiles<-trackfiles[order(trackfiles,decreasing=T)]
	
	header<-c(
		"",
		"",
		paste("track ",parentname,"-parent",sep=""),
		"compositeTrack on",
		paste("shortLabel",parentname),
		paste("longLabel",parentname),
		"type bed 3",
		"noInherit on",
		"configurable on",
		"visilbility squish",
		paste("subGroup1 view Views ",parentname,"=",parentname,sep=""),
		"",
		paste("track",parentname),
		paste("shortLabel",parentname),
		paste("longLabel",parentname),
		paste("view ",parentname),
		paste("parent ",parentname,"-parent",sep=""),
		"graphType bar",
		"configurable on",
		"visilbility squish",
		"maxHeightPixels 200:64:32",
		"autoScale off",
		"windowingFunction mean",
		"smoothingWindow off",
		""
		)

	for(i in 1:length(trackfiles)){
		c<-paste(col2rgb(plotcolors[i]),collapse=",")
		track<-c(
			paste("track ",parentname,"-",tracknames[i],sep=""),
			paste("bigDataUrl bbi/",trackfiles[i],sep=""),
			#paste("shortLabel ",letters[i],"_",tracknames[i],sep=""),
			paste("shortLabel",tracknames[i]),
			paste("longLabel",tracknames[i]),
			paste("type bigWig",range[1],range[2]),
			paste("parent",parentname),
			paste("color",c),
			paste("altColor",c),
			paste("subGroups view=",parentname,sep=""),
			"visibility squish",
			""
		)
		header<-append(header,track)
	}
	
	outfile<-data.frame("V1"=header,stringsAsFactors=FALSE)
	write.tsv(outfile,file="tmphubdb.txt")
	filelist<-paste(trackfiles,sep=" ",collapse=" ")
	login<-unlist(strsplit(hubloc,":"))[1]
	path<-unlist(strsplit(hubloc,":"))[2]
	cat(paste("scp ",hubloc,"/bbi/",sep=""))
	system(paste("scp ",filelist," ",hubloc,"/bbi/",sep=""))
	system(paste("cat tmphubdb.txt | ssh ",login," 'cat >> ",path,"/hubDb.txt'",sep=""))
}
