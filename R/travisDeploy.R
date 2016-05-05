
travisDeploy<-function(message="no message"){
	library(travis)
	res <- system(paste0("git -C /lustre/maize/home/dlv04c/software/r/travis/ add /lustre/maize/home/dlv04c/software/r/travis/ &&\
	git -C /lustre/maize/home/dlv04c/software/r/travis/ commit -a -m '",message,"' &&\
	git -C /lustre/maize/home/dlv04c/software/r/travis/ push"))
	library(devtools)
	detach("package:travis",unload=T)
	install_github("dvera/travis")
	library(travis)
}
