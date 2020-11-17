## Extract results of interest, write TAF output tables
library(icesTAF)
taf.library(FLFleet)
taf.library(FLCore)

mkdir("output")

# load project mapping
project_mapping <- read.csv("bootstrap/data/project_mapping.csv")

load(file.path("model",'HERAS_baseline.RData'))

# write NSAS index at age
write.taf(flr2taf(HERAS_baseline$HERAS.NSAS$HERAS_all@index), 
          file = "NSAS_indexAtAge.csv",
          dir="output")
# write NSAS SSB
write.taf(flr2taf(HERAS_baseline$stk.NSAS@stock), 
          file = "NSAS_SSB.csv",
          dir="output")
# write NSAS proportion mature
write.taf(flr2taf(HERAS_baseline$stk.NSAS@mat), 
          file = "NSAS_mat.csv",
          dir="output")
# write NSAS weight at age
write.taf(flr2taf(HERAS_baseline$stk.NSAS@stock.wt), 
          file = "NSAS_wa.csv",
          dir="output")



# write WBSS index at age
write.taf(flr2taf(HERAS_baseline$HERAS.WBSS$HERAS_all@index), 
          file = "WBSS_indexAtAge.csv",
          dir="output")
# write WBSS SSB
write.taf(flr2taf(HERAS_baseline$stk.WBSS@stock), 
          file = "WBSS_SSB.csv",
          dir="output")
# write WBSS proportion mature
write.taf(flr2taf(HERAS_baseline$stk.WBSS@mat), 
          file = "WBSS_mat.csv",
          dir="output")
# write WBSS weight at age
write.taf(flr2taf(HERAS_baseline$stk.NSAS@stock.wt), 
          file = "WBSS_wa.csv",
          dir="output")


# some example of StoX output (superIndividual table)
for(idxProject in 1:dim(project_mapping)[1]){
  msg("creating outout for: ", project_mapping$project_name[idxProject])
  load(file.path("model", paste0('BaseLineReport_',
                                 project_mapping$project_name[idxProject], ".RData")))
  write.taf(
    as.data.frame(currentBaseLineReport$outputData$FillMissingData),
    file = paste0(idxProject, "FillMissingData.csv"),
    dir = "output")
}
