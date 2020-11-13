fill_missing_species_v3 <- function(endTab,mode='Cells'){

  #endTab <- endTab2
  endTab$filling_weight_length <- 0
  
  #######################################################
  # find missing cells for obvious case using full table
  #######################################################
  
  # weight/length relationship
  growth <- nls(formula = as.numeric(weight) ~ a*as.numeric(length)^b,
                data = endTab[endTab$weight != '-' | endTab$length == '-', ],
                start = list(a = 0.01, b= 2.8))
  
  idxFillMat <- which(  endTab$specialstage == '-'| 
                          endTab$stage == '-' | 
                          endTab$age == '-' | 
                          endTab$weight == '-' |
                          endTab$length == '-')
  
  endTabOut <- endTab[idxFillMat,]
  
  # fill in age 0 to immature
  idxEmptyMat     <-  which(endTabOut$specialstage=="-")
  idxTemp         <- idxEmptyMat[endTabOut$age[idxEmptyMat] == 0]
  if(length(idxTemp != 0)){
    print(paste0('fill in age 0 to immature - N=',length(idxTemp)))
    endTabOut[idxTemp,]$specialstage  <- 'IMM'
    endTabOut[idxTemp,]$filling_weight_length       <- 1
  }
  
  # fill in fish of length < 8.5 cm as IMM
  idxEmptyMat   <-  which(endTabOut$specialstage=="-")
  idxTemp       <- idxEmptyMat[endTabOut$length[idxEmptyMat] < 8.5]
  if(length(idxTemp != 0)){
    print(paste0('fill in fish of length < 8.5 cm as IMM - N=',length(idxTemp)))
    endTabOut[idxTemp,]$specialstage <- 'IMM'
    endTabOut[idxTemp,]$filling_weight_length       <- 1
  }

  # fill in missing lengths
  idxEmptyLength  <-  which(endTabOut$length=="-")
  if(length(idxEmptyLength) != 0){
    print(paste0('fill in missing lengths - N=',length(idxEmptyLength)))
    endTabOut$length[idxEmptyLength] <- (endTabOut$weight[idxEmptyLength]/coef(growth)[1])^(1-coef(growth)[2])
    endTabOut$filling_weight_length[idxEmptyLength] <- 1
  }

  
  # fill in missing weights
  idxEmptyWeight  <-  which(endTabOut$weight=="-")
  if(length(idxEmptyWeight) != 0){
    print(paste0('fill in missing weights - N=',length(idxEmptyWeight)))
    yPred <- predict(object=growth,newdata=data.frame(length=endTabOut$length[idxEmptyWeight]))
    endTabOut$weight[idxEmptyWeight] <- yPred
    endTabOut$filling_weight_length[idxEmptyWeight] <- 1
  }
  
  # fill in missing ages
  idxEmptyAge  <-  which(endTabOut$age=="-")
  
  uniqueLength <- unique(endTabOut$length[idxEmptyAge])
  uniqueStrata <- unique(endTabOut$Stratum[idxEmptyAge])
  
  if(length(idxEmptyAge) != 0){
    print(paste0('fill in missing ages - N=',length(idxEmptyAge)))
    
    for(idxLength in uniqueLength){
      for(idxStrata in uniqueStrata){
        endTabFiltAge <- endTab[  endTab$Stratum ==  idxStrata&
                                    endTab$length >= idxLength-2 &
                                    endTab$length <= idxLength+2 & 
                                    endTab$age != '-',]
        endTabFiltAge <- aggregate(endTabFiltAge$Abundance,list(endTabFiltAge$age),sum)
        names(endTabFiltAge) <- c('age_group','Abundance')
        endTabOut[endTabOut$Stratum == idxStrata & endTabOut$length == idxLength,]$age <- as.numeric(endTabFiltAge[which(endTabFiltAge$Abundance == max(endTabFiltAge$Abundance)),]$age_group)
      }
    }
  }
  
  # replace values in tab
  endTab[idxFillMat,] <- endTabOut

  ###########################################
  # aggregate table
  ###########################################
  
  # calculate biomass for each entry
  endTab$biomass<-as.numeric(as.character(endTab$Abundance))*as.numeric(as.character(endTab$weight))
  
  abun <- aggregate(endTab$Abundance,
                   list(  endTab$Stratum,
                          endTab$length,
                          endTab$weight,
                          endTab$age,
                          endTab$specialstage,
                          endTab$stage,
                          endTab$filling_weight_length),
                   sum)
  
  names(abun)<-c("strata","length","weight","age","maturity","stock","filling_weight_length","Abundance")
  
  biom <- aggregate(endTab$biomass, 
                   list(  endTab$Stratum, 
                          endTab$length,
                          endTab$weight,
                          endTab$age, 
                          endTab$specialstage, 
                          endTab$stage,
                          endTab$filling_weight_length), 
                   sum)
  names(biom)<-c("strata","length","weight","age","maturity","stock","filling_weight_length","biomass_g")
  
  all<-abun
  all$biomass_g<- biom$biomass_g
  all$meanW_g<-all$biomass_g/all$Abundance
  all$filling_mat_ID <- 0
  
  endTab <- all
  
  ###########################################
  # fill in maturity and stock id
  ###########################################
  
  idxFillMat <- which(  endTab$maturity == '-'| 
                          endTab$stock == '-' | 
                          endTab$age == '-' | 
                          endTab$length == '-')
  
  endTabOut <- endTab[idxFillMat,]
  endTabOut$filling_mat_ID <- 1
  
  for(idxSpecimen in 1:dim(endTabOut)[1]){
    #print(idxSpecimen)
    missingFields <- c(endTabOut[idxSpecimen,]$stock == '-', endTabOut[idxSpecimen,]$maturity == '-')
    
    
    # if both fields missing
    if(all(missingFields == c(TRUE, TRUE))){
      
      # 1: check strata. NSAS if part of a strata other than those where mixing is expected
      if(is.na(any(match(c(151,152,41,42,31,21),endTabOut$strata[idxSpecimen])))){
        endTabOut[idxSpecimen,]$stock <- 'her-47d3'
      }else{# 2: if in strata with NSAS/WBSS mixing, assign to highest probability based on maturity and age. If not age check -+2 cm of length
        length_inter  <- 0
        probAllNSAS   <- -1
        probAllIMM    <- -1
        flag_age      <- FALSE

        # if there is an age
        if(any(endTab$strata == endTabOut$strata[idxSpecimen]&
               endTab$age == endTabOut$age[idxSpecimen]&
               endTab$stock != '-' &
               endTab$maturity != '-')){
          endTabFiltAge <- endTab[0,0]
          while(((probAllNSAS == -1 & probAllIMM == -1) | (probAllNSAS == 0.5 & probAllIMM == 0.5)) & 
                dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                             endTab$age == endTabOut$age[idxSpecimen]&
                             endTab$stock != '-' &
                             endTab$maturity != '-',])[1] > dim(endTabFiltAge)[1]){
            # create table
            endTabFiltAge <- endTab[  endTab$strata == endTabOut$strata[idxSpecimen] &
                                        endTab$age == endTabOut$age[idxSpecimen] &
                                        endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                        endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                        endTab$stock != '-' &
                                        endTab$maturity != '-',]
            length_inter <- length_inter + 1
            
            if(dim(endTabFiltAge)[1] != 0){
              # calculate probabilities
              if(mode == 'Cells'){
                probAllNSAS <- dim(endTabFiltAge[endTabFiltAge$stock == 'her-47d3',])[1]/dim(endTabFiltAge)[1]
                probAllIMM  <- dim(endTabFiltAge[endTabFiltAge$maturity == 'IMM',])[1]/dim(endTabFiltAge)[1]
              }else if(mode == 'Abundance'){
                probAllNSAS <- sum(endTabFiltAge[endTabFiltAge$stock == 'her-47d3',]$Abundance)/sum(endTabFiltAge$Abundance)
                probAllIMM  <- sum(endTabFiltAge[endTabFiltAge$maturity == 'IMM',]$Abundance)/sum(endTabFiltAge$Abundance)
              }
            }
          }
          if(probAllNSAS != 0.5 & probAllIMM != 0.5)flag_age <- TRUE
        }

        # if no age, go with length
        if(flag_age == FALSE){
          endTabFiltLength <- endTab[0,0]
          while(((probAllNSAS == -1 & probAllIMM == -1) | (probAllNSAS == 0.5 & probAllIMM == 0.5)) &
                dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                             endTab$age == endTabOut$age[idxSpecimen]&
                             endTab$stock != '-' &
                             endTab$maturity != '-',])[1] > dim(endTabFiltLength)[1]){
            # create table
            endTabFiltLength <- endTab[ endTab$strata == endTabOut$strata[idxSpecimen] &
                                          endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                          endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                          endTab$stock != '-' &
                                          endTab$maturity != '-',]
            
            # calculate probabilities
            if(dim(endTabFiltLength)[1] != 0){
              if(mode == 'Cells'){
                probAllNSAS <- dim(endTabFiltLength[endTabFiltLength$stock == 'her-47d3',])[1]/dim(endTabFiltLength)[1]
                probAllIMM  <- dim(endTabFiltLength[endTabFiltLength$maturity == 'IMM',])[1]/dim(endTabFiltLength)[1]
              }else if(mode == 'Abundance'){
                probAllNSAS <- sum(endTabFiltLength[endTabFiltLength$stock == 'her-47d3',]$Abundance)/sum(endTabFiltLength$Abundance)
                probAllIMM  <- sum(endTabFiltLength[endTabFiltLength$maturity == 'IMM',]$Abundance)/sum(endTabFiltLength$Abundance)
              }
              length_inter <- length_inter + 1
            }else length_inter <- length_inter + 1
          }
        }
        
        if(probAllNSAS >= 0.5)endTabOut[idxSpecimen,]$stock <- 'her-47d3'else endTabOut[idxSpecimen,]$stock <- 'her-3a22'
        if(probAllIMM >= 0.5)endTabOut[idxSpecimen,]$maturity <- 'IMM'else endTabOut[idxSpecimen,]$maturity <- 'MAT'
        
      }
      
      
    # if missing stock ID only
    }else if(which(missingFields) == 1){
      # 1: check strata. NSAS if part of a strata other than those where mixing is expected
      if(is.na(any(match(c(151,152,41,42,31,21),endTabOut$strata[idxSpecimen])))){
        endTabOut[idxSpecimen,]$stock <- 'her-47d3'
        
      }else{# 2: if in strata with NSAS/WBSS mixing, assign to highest probability based on maturity and age. If not age check -+2 cm of length
        length_inter  <- 0
        probAllNSAS   <- -1
        probAllIMM    <- -1
        flag_age      <- FALSE

        # if there is an age
        if(any(endTab$strata == endTabOut$strata[idxSpecimen]&
               endTab$maturity == endTabOut$maturity[idxSpecimen] &
               endTab$age == endTabOut$age[idxSpecimen]&
               endTab$stock != '-')){
          endTabFiltAge <- endTab[0,0]
          while((probAllNSAS == -1 | probAllNSAS == 0.5) & 
                dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                             endTab$maturity == endTabOut$maturity[idxSpecimen] &
                             endTab$age == endTabOut$age[idxSpecimen]&
                             endTab$stock != '-',])[1] > dim(endTabFiltAge)[1]){
            endTabFiltAge <- endTab[  endTab$strata == endTabOut$strata[idxSpecimen] &
                                        endTab$maturity == endTabOut$maturity[idxSpecimen] &
                                        endTab$age == endTabOut$age[idxSpecimen] &
                                        endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                        endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                        endTab$stock != '-',]
            length_inter <- length_inter + 1
            
            if(dim(endTabFiltAge)[1] != 0){
              if(mode == 'Cells'){
                probAllNSAS <- dim(endTabFiltAge[endTabFiltAge$stock == 'her-47d3',])[1]/dim(endTabFiltAge)[1]
              }else if(mode == 'Abundance'){
                probAllNSAS <- sum(endTabFiltAge[endTabFiltAge$stock == 'her-47d3',]$Abundance)/sum(endTabFiltAge$Abundance)
              }
            }
          }
          if(probAllNSAS != 0.5)flag_age <- TRUE
        }
        
        # if no age, go with length
        if(flag_age == FALSE){
          endTabFiltLength <- endTab[0,0]
          while((probAllNSAS == -1 | probAllNSAS == 0.5) &
                dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                             endTab$maturity == endTabOut$maturity[idxSpecimen] &
                             endTab$stock != '-',])[1] > dim(endTabFiltLength)[1]){
            endTabFiltLength <- endTab[ endTab$strata == endTabOut$strata[idxSpecimen] &
                                          endTab$maturity == endTabOut$maturity[idxSpecimen] &
                                          endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                          endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                          endTab$stock != '-',]
            if(dim(endTabFiltLength)[1] != 0){
              if(mode == 'Cells'){
                probAllNSAS <- dim(endTabFiltLength[endTabFiltLength$stock == 'her-47d3',])[1]/dim(endTabFiltLength)[1]
              }else if(mode == 'Abundance'){
                probAllNSAS <- sum(endTabFiltLength[endTabFiltLength$stock == 'her-47d3',]$Abundance)/sum(endTabFiltLength$Abundance)
              }
              length_inter <- length_inter + 1
            }else length_inter <- length_inter + 1
          }
        }

        # allocate species
        if(probAllNSAS >= 0.5)endTabOut[idxSpecimen,]$stock <- 'her-47d3'else endTabOut[idxSpecimen,]$stock <- 'her-3a22'
      }
      
      # if missing stock maturity only
    }else if(which(missingFields) == 2){
      length_inter  <- 0
      probAllNSAS   <- -1
      probAllIMM    <- -1
      flag_age      <- FALSE

      # if there is an age
      if(any(endTab$strata == endTabOut$strata[idxSpecimen]&
             endTab$stock == endTabOut$stock[idxSpecimen] &
             endTab$age == endTabOut$age[idxSpecimen]&
             endTab$maturity != '-')){
        endTabFiltAge <- endTab[0,0]
        while((probAllIMM == -1 | probAllIMM == 0.5) & 
              dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                           endTab$stock == endTabOut$stock[idxSpecimen] &
                           endTab$age == endTabOut$age[idxSpecimen]&
                           endTab$maturity != '-',])[1] > dim(endTabFiltAge)[1]){
          endTabFiltAge <- endTab[  endTab$strata == endTabOut$strata[idxSpecimen] &
                                      endTab$stock == endTabOut$stock[idxSpecimen] &
                                      endTab$age == endTabOut$age[idxSpecimen] &
                                      endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                      endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                      endTab$maturity != '-',]
          length_inter <- length_inter + 1
          
          if(dim(endTabFiltAge)[1] != 0){
            if(mode == 'Cells'){
              probAllIMM  <- dim(endTabFiltAge[endTabFiltAge$maturity == 'IMM',])[1]/dim(endTabFiltAge)[1]
            }else if(mode == 'Abundance'){
              probAllIMM  <- sum(endTabFiltAge[endTabFiltAge$maturity == 'IMM',]$Abundance)/sum(endTabFiltAge$Abundance)
            }
          }
        }
        if(probAllIMM != 0.5)flag_age <- TRUE
      }

      # if no age, go with length
      if(flag_age == FALSE){
        endTabFiltLength <- endTab[0,0]
        while((probAllIMM == -1 | probAllIMM == 0.5)&
              dim(endTab[  endTab$strata == endTabOut$strata[idxSpecimen]&
                           endTab$stock == endTabOut$maturity[idxSpecimen] &
                           endTab$maturity != '-',])[1] > dim(endTabFiltLength)[1]){
          endTabFiltLength <- endTab[ endTab$strata == endTabOut$strata[idxSpecimen] &
                                        endTab$stock == endTabOut$stock[idxSpecimen] &
                                        endTab$length >= (endTabOut$length[idxSpecimen]-length_inter) &
                                        endTab$length <= (endTabOut$length[idxSpecimen]+length_inter) &
                                        endTab$maturity != '-',]
          if(dim(endTabFiltLength)[1] != 0){
            if(mode == 'Cells'){
              probAllIMM  <- dim(endTabFiltLength[endTabFiltLength$maturity == 'IMM',])[1]/dim(endTabFiltLength)[1]
            }else if(mode == 'Abundance'){
              probAllIMM  <- sum(endTabFiltLength[endTabFiltLength$maturity == 'IMM',]$Abundance)/sum(endTabFiltLength$Abundance)
            }
            length_inter <- length_inter + 1
          }else length_inter <- length_inter + 1
        }
      }
      
      # allocate maturity field
      if(probAllIMM >= 0.5)endTabOut[idxSpecimen,]$maturity <- 'IMM'else endTabOut[idxSpecimen,]$maturity <- 'MAT'
    }
  }
  
  endTab[idxFillMat,] <- endTabOut
  
  return(endTab)
}