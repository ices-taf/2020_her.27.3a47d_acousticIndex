## Preprocess data, write TAF data tables

library(icesTAF)

# create data directory
mkdir("data")


data_mapping <- read.csv("bootstrap/data/data_mapping.csv")

# creat stox folder tree

project_years <- unique(data_mapping$year)

for(idxYear in project_years){
  currentProjects <- unique(data_mapping[data_mapping$year == idxYear,]$project)

  for(idxProject in currentProjects){
    wkdir <- file.path("data", "stox-project", as.character(idxYear), as.character(idxProject))
    mkdir(wkdir)

    # copy polygon files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'polygons',]$file)
    file.copy(from=file.path("bootstrap", "data", filestocopy),
              to=wkdir,
              copy.mode = TRUE)

    # create input folder
    mkdir(file.path(wkdir, 'input', 'acoustic'))
    mkdir(file.path(wkdir, 'input', 'biotic'))

    # copy acoustic files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'acoustic',]$file)
    for(idxFile in filestocopy){
      file.copy(from=file.path('bootstrap', 'data', idxFile),
                to=file.path(wkdir, "input", 'acoustic'),
                copy.mode = TRUE)
    }

    # copy biotic files
    filestocopy <- as.character(data_mapping[ data_mapping$year == idxYear &
                                                data_mapping$project == idxProject &
                                                data_mapping$type == 'biotic',]$file)
    for(idxFile in filestocopy){
      file.copy(from=file.path('bootstrap', 'data', idxFile),
                to=file.path(wkdir, 'input', 'biotic'),
                copy.mode = TRUE)
    }

    # create model and output folder
    mkdir(file.path(wkdir, 'output'))
    mkdir(file.path(wkdir, 'process'))

    # copy project file and rename it
    filestocopy <- paste0(as.character(idxProject),'.xml')
    file.copy(from=file.path('bootstrap', 'data', filestocopy),
              to=file.path(wkdir, 'process', 'project.xml'),
              copy.mode = TRUE)
  }
}
