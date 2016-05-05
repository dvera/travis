mat.plotsamplerows <-
function( matname,samples=5,ylims=c(0,30)){
	mat<-read.mat(matname)
	rows<-sample(1:nrow(mat),samples)
	plotcolors=rainbow(samples)
	plot(1:ncol(mat),mat[rows[1],],col=plotcolors[1],type="l",ylim=ylims)
	for(i in 2:samples){
		lines(1:ncol(mat),mat[rows[i],],col=plotcolors[i])
	}
	legend("topright",legend=row.names(mat)[rows],col=plotcolors,lwd=2)
}
