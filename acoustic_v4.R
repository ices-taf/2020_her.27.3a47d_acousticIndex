
library(httr)
library(jsonlite)

# These are Acoustic WebAPI for downloading StoX files:
# 1.	List of all submission -> http://acoustic.ices.dk/services/Download/GetUploadlist/
# This service returns list of submitted cruise. You can get id for each cruise and pass to the bellow services to download Acoustic or Biotic StoX files.
#
# 2.	Download Acoustic StoX file -> http://acoustic.ices.dk/services/Download/DownloadStoXAcousticXML/{ID}
# eg. http://acoustic.ices.dk/services/Download/DownloadStoXAcousticXML/3
#
# 3.	Download Biotic StoX file ->  http://acoustic.ices.dk/services/Download/DownloadStoXBioticXML/{ID}
# eg. http://acoustic.ices.dk/services/Download/DownloadStoXBioticXML/3

# -------------------
# get list of uploads
# -------------------

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
    id = uploads$cruiseField.id,
    survey = unname(unlist(uploads$cruiseField.surveyField)),
    startDate = uploads$cruiseField.startDateField,
    workingGroup = uploads$cruiseField.workingGroup,
    country = uploads$cruiseField.countryField.iDREFField,
    platform = uploads$cruiseField.platformField.iDREFField
  )

head(dat)

# -------------------
# get biota files
# -------------------

GET_biotic <- function(id, what = c("Biotic", "Acoustic"), ...) {
  what <- match.arg(what)
  res <-
    GET(
      sprintf(
        "https://acoustic.ices.dk/Submissions/Download%sCSV/%s",
        what,
        id
      )
    )
  
  tmpfile <- tempfile(fileext = ".zip")
  writeBin(content(res, as = "raw"), tmpfile)
  files <- unzip(tmpfile, list = TRUE)
  filetoget <- grep("*.csv", files$Name, value = TRUE)
  unzip(
    tmpfile,
    files = filetoget,
    exdir = tempdir()
  )
  
  unlink(tmpfile)
  
  # read in
  vars <- c("Cruise", "Haul", "Catch", "Biology")
  rows <- as.vector(readLines(file.path(tempdir(), filetoget)))
  data <- data.frame(col1=rows)
  data <- data.frame(do.call('rbind', strsplit(as.character(data$col1),',',fixed=TRUE)))
  tab <- as.character(unique(data$X1))
  
  out <- list()
  for(j in tab){
    tmp <- subset(data,X1==j)
    colnames(tmp) <- as.vector(as.character(t(tmp[1,])))
    tmp <- tmp[-1,]
    assign(j,tmp)
    rm(tmp)
    out[[which(tab==j)]] <- get(j)
  }
  
  
  return(out)
}

ids <- dat$id[dat$survey == "AC_Survey_HERAS"]

tmp <- GET_biotic(ids[1]) ### get all data frames in one list
head(tmp)





#biotics <- list()
#for (i in ids[1:2]) {
#  biotics <- c(biotics, i = list(GET_biotic(i)))
#}



## process files
