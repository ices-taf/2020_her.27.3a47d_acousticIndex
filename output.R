## Extract results of interest, write TAF output tables
library(icesTAF)

mkdir("output")

data_mapping <- read.csv("bootstrap/data/data_mapping.csv")
project_years <- unique(data_mapping$year)

for(idxYear in project_years){
  currentProjects <- as.character(unique(data_mapping[data_mapping$year == idxYear,]$project))

  for(idxProject in currentProjects){
    msg("creating outout for: ", idxProject)
    load(file.path("model", paste0(idxProject, "_BaseLineReport.RData")))

    # test output creation
    write.taf(
      as.data.frame(currentBaseLineReport$outputData$FillMissingData),
      file = paste0(idxProject, "FillMissingData.csv"),
      dir = "output"
    )
  }
}
