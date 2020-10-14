## Run analysis, write model results

# run these set up options first
options(java.parameters = "-Xmx2000m") # increase RAM for Java VM
library(rJava)
.jinit()

library(icesTAF)
taf.library(Rstox)
library(mgcv)

# create model dir
mkdir("model")

# copy over data (create stucture for stox in model folder)
cp(from = 'data/stox-project', to = 'model')


data_mapping <- read.csv("bootstrap/data/data_mapping.csv")
project_years <- unique(data_mapping$year)

# get configuration file
for (idxYear in project_years) {
  currentProjects <- as.character(unique(data_mapping[data_mapping$year == idxYear, ]$project))

  for (idxProject in currentProjects) {
    msg("running Baseline for: ", idxProject)

    wkdir <- file.path(getwd(), "model", "stox-project", as.character(idxYear), as.character(idxProject))

    openProject(wkdir)

    runBaseline(projectName = wkdir, save = TRUE, exportCSV = TRUE, modelType = c("baseline"))
    runBaseline(projectName = wkdir, save = TRUE, exportCSV = TRUE, modelType = c("baseline-report"))

    currentBaseLineReport <- getBaseline(wkdir, save = FALSE, exportCSV = FALSE, modelType = c("baseline-report"))

    save(currentBaseLineReport,
      file = file.path("model", paste0(idxProject, "_BaseLineReport.RData")),
      compress = "xz"
    )

    # StoX bootstraping
    # runBootstrap(projectName=file.path(getwd(),idxProject),
    #             bootstrapMethod="AcousticTrawl",
    #             acousticMethod="PSU~Stratum",
    #             bioticMethod="PSU~Stratum",
    #             startProcess="TotalLengthDist",
    #             endProcess="SuperIndAbundance",
    #             nboot=nboot,
    #             seed=1, cores=4)

    # saveProjectData(dataDirs[idxDataDir])

    # garbage collect
    .jcall("java/lang/System", method = "gc")
  }
}
