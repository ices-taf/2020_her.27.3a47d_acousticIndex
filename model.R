## Run analysis, write model results

# run these set up options first
options(java.parameters = "-Xmx2000m") # increase RAM for Java VM
library(rJava)
.jinit()

library(icesTAF)
taf.library(Rstox)
taf.library(FLFleet)
taf.library(FLCore)
library(mgcv)

# source functions
source('bootstrap/software/utilities/fill_missing_mat_v1.R')
source('bootstrap/software/utilities/compute_NSAS_WBSS_strata.R')
source('bootstrap/software/utilities/fill_missing_species_v3.R')

# create model dir
mkdir("model")

# parameters
surveyYear <- 2019
minYear <- surveyYear
maxYear <- surveyYear
yearVec <- minYear:maxYear
nYears  <- length(yearVec)
minAge  <- 0
maxAge  <- 9
ageVec  <- minAge:maxAge
nAges   <- length(ageVec)
nIter   <- 1


# load NO split
load('data/split_prop_NO.Rdata')

# copy over data (create stucture for stox in model folder)
cp(from = 'data/stox-project', to = 'model')

# load project mapping
project_mapping <- read.csv("bootstrap/data/project_mapping.csv")

######################################################################
######################################################################
######################################################################
#### # initialise FLStock objects
######################################################################
######################################################################
######################################################################

# stock object
FLR_temp <- FLQuant( array(-1,dim = c(nAges,nYears,1,1,1,nIter)),
                     dimnames=list( age=as.character(ageVec),
                                    year=as.character(yearVec),
                                    unit=c("1e6"),
                                    season="all", 
                                    area='all'),
                     units='1e6')

NSAS <- FLStock(name = 'NSAS', 
                desc = 'NSAS',
                FLR_temp)

WBSS <- FLStock(name = 'WBSS', 
                desc = 'WBSS',
                FLR_temp)

######################################################################
#### # initialise FLIndex objects
######################################################################

# EU NO index object
FLR_HERAS_EU_NO <- FLQuant( array(-1,dim = c(nAges,nYears,1,1,2,nIter)),
                            dimnames=list(age=as.character(ageVec),
                                          year=as.character(yearVec),
                                          unit=c("1e6"),
                                          season="all", 
                                          area=c('EU','NO')),
                            units='1e6')
FLR_HERAS_EU_NO <- FLIndex(name = 'HERAS_EU_NO', 
                           desc = 'HERAS index from EU and NO StoX projects', 
                           index=FLR_HERAS_EU_NO)

# total index object
FLR_HERAS <- FLQuant( array(-1,dim = c(nAges,nYears,1,1,1,nIter)),
                      dimnames=list(age=as.character(ageVec),
                                    year=as.character(yearVec),
                                    unit="1e6",
                                    season="all", 
                                    area='all'),
                      units='1e6')
FLR_HERAS <- FLIndex( name = 'HERAS_all', 
                      desc = 'HERAS index (total)', 
                      index=FLR_HERAS)

HERAS.NSAS <- FLIndices(FLR_HERAS_EU_NO, FLR_HERAS)
HERAS.WBSS <- FLIndices(FLR_HERAS_EU_NO, FLR_HERAS)

######################################################################
# run split projects
######################################################################

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
  #             nboot=nIter,
  #             seed=1, cores=4)
  
  #saveProjectData(wkdir)
  
  .jcall("java/lang/System", method = "gc")
}

#################################################################
## create index and assign values
#################################################################
outIndex <- compute_NSAS_WBSS_strata(   endTabEU,
                                        endTabNO,
                                        split_prop,
                                        ageVec,
                                        surveyYear)

# assign numbers at age to FLR objects
HERAS.NSAS$HERAS_EU_NO@index[,ac(surveyYear),,,'NO'] <- colSums(outIndex$NSAS_NO)*1e-6
HERAS.WBSS$HERAS_EU_NO@index[,ac(surveyYear),,,'NO'] <- colSums(outIndex$WBSS_NO)*1e-6
HERAS.NSAS$HERAS_EU_NO@index[,ac(surveyYear),,,'EU'] <- colSums(outIndex$NSAS_EU)*1e-6
HERAS.WBSS$HERAS_EU_NO@index[,ac(surveyYear),,,'EU'] <- colSums(outIndex$WBSS_EU)*1e-6

# fill in stock weight at age and maturity
NSAS@stock.wt[,ac(surveyYear)]  <- outIndex$wa$NSAS
NSAS@mat[,ac(surveyYear)]       <- outIndex$propM$NSAS
NSAS@stock.n[,ac(surveyYear)]   <- areaSums(HERAS.NSAS$HERAS_EU_NO@index[,ac(surveyYear)])

WBSS@stock.wt[,ac(surveyYear)]  <- outIndex$wa$WBSS
WBSS@mat[,ac(surveyYear)]       <- outIndex$propM$WBSS
WBSS@stock.n[,ac(surveyYear)]   <- areaSums(HERAS.WBSS$HERAS_EU_NO@index[,ac(surveyYear)])

# fill in SSB
NSAS@stock[,ac(surveyYear)] <- sum(NSAS@stock.n[,ac(surveyYear)]*NSAS@stock.wt[,ac(surveyYear)]*NSAS@mat[,ac(surveyYear)])

WBSS@stock[,ac(surveyYear)] <- sum(WBSS@stock.n[,ac(surveyYear)]*WBSS@stock.wt[,ac(surveyYear)]*WBSS@mat[,ac(surveyYear)])

# combine NO and EU components
HERAS.NSAS$HERAS_all@index[] <- areaSums(HERAS.NSAS$HERAS_EU_NO@index)
HERAS.WBSS$HERAS_all@index[] <- areaSums(HERAS.WBSS$HERAS_EU_NO@index)

HERAS_baseline <- list(HERAS.NSAS=HERAS.NSAS,
                       HERAS.WBSS=HERAS.WBSS,
                       stk.WBSS=WBSS,
                       stk.NSAS=NSAS)

save(HERAS_baseline,
     file = file.path("model",'HERAS_baseline.RData'))
