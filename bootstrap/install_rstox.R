library(icesTAF)

mkdir("../../library")

dep.pck <- c(	"data.table", 
				"ggplot2", 
				"pbapply", 
				"rgdal", 
				"rgeos", 
				"rJava", 
				"sp", 
				"XML",
				"mgcv")
install.packages(dep.pck, repos="http://cran.us.r-project.org", type="binary",lib="../../library")

install.packages("ftp://ftp.imr.no/StoX/Download/Rstox/Rstox_1.11.tar.gz", repos=NULL,  lib="../../library")
