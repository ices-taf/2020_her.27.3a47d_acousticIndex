## Run analysis, write model results

# run these set up options first
options(java.parameters = "-Xmx2000m") # increase RAM for Java VM
library(rJava)
.jinit()

library(icesTAF)
taf.library(Rstox)
library(mgcv)

# source functions
source('bootstrap/software/utilities/fill_missing_mat_v1.R')
source('bootstrap/software/utilities/compute_NSAS_WBSS_strata.R')
source('bootstrap/software/utilities/fill_missing_species_v3.R')

# create model dir
mkdir("model")

# age parameters
minAge  <- 0
maxAge  <- 9
ageVec  <- minAge:maxAge


# load NO split
load('data/split_prop_NO.Rdata')

# copy over data (create stucture for stox in model folder)
cp(from = 'data/stox-project', to = 'model')

# load project mapping
project_mapping <- read.csv("bootstrap/data/project_mapping.csv")

# run split projects
for(idxProject in 1:dim(project_mapping)[1]){
  msg("running Baseline for: ", project_mapping$project_name[idxProject])
  
  wkdir <- file.path(getwd(), "model", "stox-project", as.character(project_mapping$project_name[idxProject]))
  
  openProject(wkdir)
  runBaseline(projectName = wkdir, save = TRUE, exportCSV = TRUE, modelType = c("baseline"))
  runBaseline(projectName = wkdir, save = TRUE, exportCSV = TRUE, modelType = c("baseline-report"))
  
  currentBaseLineReport <- getBaseline(wkdir, save = FALSE, exportCSV = FALSE, modelType = c("baseline-report"))
  
  save(currentBaseLineReport,
       file = file.path("model", paste0('BaseLineReport_',
                                        project_mapping$project_name[idxProject], ".RData")),
       compress = "xz")
  
  if(project_mapping$component[idxProject] == 'EU'){
    endTabEU <- currentBaseLineReport$outputData$FillMissingData
    endTabEU[is.na(endTabEU)] <- '-' # replace NA with '-'
  }else if((project_mapping$component[idxProject] == 'NO')){
    endTabNO <- currentBaseLineReport$outputData$FillMissingData
    endTabNO[is.na(endTabNO)] <- '-' # replace NA with '-'
  }
  
  # bootstraping
  #runBootstrap(projectName=wkdir,
  #             bootstrapMethod="AcousticTrawl",
  #             acousticMethod="PSU~Stratum",
  #             bioticMethod="PSU~Stratum",
  #             startProcess="TotalLengthDist",
  #             endProcess="SuperIndAbundance",
  #             nboot=500,
  #             seed=1, cores=4)
  
  #saveProjectData(wkdir)
  
  .jcall("java/lang/System", method = "gc")
}

#################################################################
## create index and assign values
#################################################################
#outIndex <- compute_NSAS_WBSS_strata(   endTabEU,
#                                        endTabNO,
#                                        split_prop,
#                                        ageVec,
#                                        surveyYear)
