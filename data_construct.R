# Preprocess data, write TAF data tables

library(icesTAF)
library(httr)
library(jsonlite)
library(dplyr)

# create data directory
mkdir("data")
mkdir("data/stox-project")

# get ICES DB uploads
res <-
  GET(
    "https://acoustic.ices.dk/services/Download/GetUploadlist/",
    accept_json()
  )
uploads <-
  parse_json(
    content(res, as = "text"),
    simplifyVector = TRUE,
    flatten = TRUE
  )

dat <-
  data.frame(
    id = uploads$Cruise.ID,
    survey = unname(unlist(uploads$Cruise.Survey)),
    startDate = uploads$Cruise.StartDate,
    workingGroup = uploads$Cruise.WorkingGroup,
    country = uploads$Cruise.Country.IDREF,
    platform = uploads$Cruise.Platform.IDREF
  )

dat$year    <- as.integer(substring(dat$startDate,1,4))
dat$country <- substring(dat$country,10,15)

# load data mapping
data_mapping <- read.csv("bootstrap/data/data_mapping.csv")

projects <- unique(data_mapping$project)

for(idxProject in projects){
  
  wkdir <- file.path("data", "stox-project", as.character(idxProject))
  mkdir(wkdir)
  
  #############################
  # creat stox folder tree
  #############################
  
  # create input folder
  mkdir(file.path(wkdir, 'input', 'acoustic'))
  mkdir(file.path(wkdir, 'input', 'biotic'))
  
  # create process and output folder
  mkdir(file.path(wkdir, 'output'))
  mkdir(file.path(wkdir, 'process'))
  
  #############################
  # copy local files
  #############################
  
  fileLocal <- data_mapping[data_mapping$project_name == idxProject & data_mapping$location == 'local',]
  
  for(idxLocal in 1:dim(fileLocal)[1]){
    if(fileLocal$type[idxLocal] == 'acoustic'){
      tempDir <- file.path(wkdir, 'input', 'acoustic')
    }else if(fileLocal$type[idxLocal] == 'biotic'){
      tempDir <- file.path(wkdir, 'input', 'biotic')
    }else if(fileLocal$type[idxLocal] == 'project'){
      tempDir <- file.path(wkdir, 'process','project.xml')
    }else if(fileLocal$type[idxLocal] == 'strata'){
      tempDir <- file.path(wkdir)
    }
    
    file.copy(from=file.path("bootstrap", "data", fileLocal$file_input[idxLocal]),
              to=tempDir,
              copy.mode = TRUE,
              overwrite = TRUE)
  }
  
  #############################
  # pull files from ICES DB
  #############################
  
  fileDB <- data_mapping[data_mapping$project_name == idxProject & data_mapping$location == 'ices_DB',]
  
  for(idxDB in 1:dim(fileDB)[1]){
    idxDB
    # set up copy folder
    if(fileDB$type[idxDB] == 'acoustic'){
      tempDir <- file.path(wkdir, 'input', 'acoustic',fileDB$file_input[idxDB])
    }else if(fileDB$type[idxDB] == 'biotic'){
      tempDir <- file.path(wkdir, 'input', 'biotic',fileDB$file_input[idxDB])
    }
    
    # find the right DB ID
    ID_DB <- dat$id[  dat$survey == "AC_Survey_HERAS" & 
                        dat$year == fileDB$year[idxDB] & 
                        dat$country == fileDB$country[idxDB]]
    
    # build request
    if(fileDB$type[idxDB] == 'acoustic'){
      res <-
        GET(  sprintf(paste0("http://acoustic.ices.dk/services/Download/DownloadStoXAcousticXML/",ID_DB),
              'acoustic',
              ID_DB))
    }else if(fileDB$type[idxDB] == 'biotic'){
      res <-
        GET(  sprintf(paste0("http://acoustic.ices.dk/services/Download/DownloadStoXBioticXML/",ID_DB),
              'biotic',
              ID_DB))
    }
    
    # pull request
    tmpfile <- tempfile(fileext = ".zip")
    tmpPath <- dirname(tmpfile)
    writeBin(content(res, as = "raw"), tmpfile)
    files <- unzip(tmpfile, list = TRUE)
    
    filetoget <- grep("*.xml", files$Name, value = TRUE)
    unzip(  tmpfile,
            files = filetoget,
            exdir = tmpPath)

    file.copy(from=file.path(tmpPath,filetoget),
              to=tempDir,
              copy.mode = TRUE,
              overwrite = TRUE)
  }
}

#############################
# Compute split for NO
#############################
load("bootstrap/data/NO_split_HERAS_2014_2020.Rdata")

HERAS <- HERAS %>% filter(lengthsamplecount>=30)
HERAS <- HERAS %>% filter(longitudestart>=2 & longitudestart<=6)

split_prop <- HERAS %>% filter(Age_Group>0 & vertebrae>0) %>% group_by(year, Stratum, Age_Group) %>% 
  summarise(Mean_VS = mean(vertebrae),
            No = length(vertebrae))

split_prop$Prop_WBSS<-(56.5-split_prop$Mean_VS)/(56.5-55.8)
split_prop$Prop_WBSS[split_prop$Prop_WBSS<0]<-0
split_prop$Prop_WBSS[split_prop$Prop_WBSS>1]<-1
split_prop$Prop_NSAS<-1-split_prop$Prop_WBSS

save(split_prop,file = 'data/split_prop_NO.Rdata')
