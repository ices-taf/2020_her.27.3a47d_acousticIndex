## Run analysis, write model results

library(icesTAF)
taf.library(Rstox)
taf.library(mgcv)
library(Rstox)

data_mapping <- read.csv("bootstrap/data/data_mapping.csv")

project_years <- unique(data_mapping$year)

nboot <- 500

setwd('stox_project')

for(idxYear in project_years){
  setwd(as.character(idxYear))
  
  currentProjects <- as.character(unique(data_mapping[data_mapping$year == idxYear,]$project))
  for(idxProject in currentProjects){
    runBaseline(projectName = file.path(getwd(),idxProject), save = TRUE,exportCSV = TRUE,modelType = c('baseline'))
  }
}