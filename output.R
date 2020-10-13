## Extract results of interest, write TAF output tables
library(icesTAF)

mkdir("output")

data_mapping <- read.csv("bootstrap/data/data_mapping.csv")

setwd('output')

for(idxYear in project_years){
	mkdir(as.character(idxYear))
  setwd(as.character(idxYear))
  
  currentProjects <- as.character(unique(data_mapping[data_mapping$year == idxYear,]$project))
  for(idxProject in currentProjects){
    load(file.path("../../model",paste0(idxProject,'_BaseLineReport.RData')))
    
    write.taf(as.data.frame(currentBaseLineReport$outputData$FillMissingData),
              "FillMissingData.csv")
  }
}

setwd('../..')