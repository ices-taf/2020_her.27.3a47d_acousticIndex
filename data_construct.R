## Preprocess data, write TAF data tables

library(icesTAF)

#rm(list=ls())

#setwd('C:/git/2020_her.27.3a47d_acousticIndex')

#library(Rstox)
#library(mgcv)

data_mapping <- read.csv("bootstrap/data/data_mapping.csv")


# creat stox folder tree



mkdir("model")

setwd('model')

project_years <- unique(data_mapping$year)

for(idxYear in project_years){
  mkdir(as.character(idxYear))
  setwd(as.character(idxYear))
  
  currentProjects <- unique(data_mapping[data_mapping$year == idxYear,]$project)
  
  for(idxProject in currentProjects){
    mkdir(as.character(idxProject))
    setwd(as.character(idxProject))
    
    # copy polygon files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'polygons',]$file)
    
    file.copy(from=paste0('../../../bootstrap/data/',filestocopy),
              to=getwd(),
              copy.mode = TRUE)
    
    # create input folder
    mkdir('input')
    setwd('input')
    mkdir('acoustic')
    mkdir('biotic')
    
    # copy acoustic files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'acoustic',]$file)
    
    for(idxFile in filestocopy){
      file.copy(from=paste0('../../../../bootstrap/data/',idxFile),
                to='acoustic',
                copy.mode = TRUE)
    }
    
    # copy biotic files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'biotic',]$file)
    
    for(idxFile in filestocopy){
      file.copy(from=paste0('../../../../bootstrap/data/',idxFile),
                to='biotic',
                copy.mode = TRUE)
    }
    
    setwd('../')
    
    
    mkdir('output')
    mkdir('process')
    
    # copy project file and rename it
    setwd('process')
    filestocopy <- paste0(as.character(idxProject),'.xml')
    file.copy(from=paste0('../../../../bootstrap/data/',filestocopy),
              to=getwd(),
              copy.mode = TRUE)
    file.rename(filestocopy,'project.xml')
    
    # going back to project root
    setwd('../../')
  }
  
  # going to root for years
  setwd('..')
}

setwd('..')