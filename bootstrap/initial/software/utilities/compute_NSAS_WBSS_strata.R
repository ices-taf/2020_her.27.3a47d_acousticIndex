compute_NSAS_WBSS_strata <- function( endTabEU,
                                      endTabNO,
                                      split_prop,
                                      ageVec,
                                      surveyYear){
  
  nAges   <- length(ageVec)
  minAge  <- min(ageVec)
  maxAge  <- max(ageVec)
  
  #################################################################
  ## NO project
  #################################################################
  
  endTabNO <- fill_missing_mat_v1(endTabNO,mode='Abundance')
  endTabNO.NSAS <- endTabNO
  endTabNO.WBSS <- endTabNO
  
  uniqueStrataNO <- unique(endTabNO$strata)
  
  # initialize array nStrata x nAges
  indexHERASComponent_NO.NSAS <- array(0, dim=c(length(uniqueStrataNO),nAges))
  indexHERASComponent_NO.WBSS <- array(0, dim=c(length(uniqueStrataNO),nAges))
  
  # split numbers at age per strata (WBSS/NSAS)
  for(idxStrata in 1:length(uniqueStrataNO)){
    strataCurrent <- uniqueStrataNO[idxStrata]
    
    # in the NO 2016 project, there is strata division (11a, 11b, 141a, 141b, 141c)
    # in the NO 2015 project, there is weird strata naming (09_NOR_South, 10_NOR_North)
    if(is.character(strataCurrent)){
      if(strataCurrent == '09_NOR_South'){
        strataCurrentNb <- 141
      }else if(strataCurrent == '10_NOR_North'){
        strataCurrentNb <- 11
      }else{
        strataCurrentNb <- as.numeric(substr(strataCurrent,1,nchar(strataCurrent)-1))
      }
    }else{
      strataCurrentNb <- strataCurrent
    }
    
    split_propArray <- as.data.frame(array(0,dim=c(nAges,3)))
    colnames(split_propArray) <- c('age','Prop_WBSS','Prop_NSAS')
    split_propArray$age <- ageVec
    
    split_propFilt  <- split_prop[split_prop$year == surveyYear & split_prop$Stratum == strataCurrentNb,]

    split_propArray[match(split_propFilt$Age_Group,split_propArray$age),]$Prop_WBSS <- split_propFilt$Prop_WBSS
    split_propArray[match(split_propFilt$Age_Group,split_propArray$age),]$Prop_NSAS <- split_propFilt$Prop_NSAS
    # all ages below minimum is set to minimum
    split_propArray[split_propArray$age < min(split_propFilt$Age_Group),]$Prop_WBSS <- 
      split_propFilt$Prop_WBSS[split_propFilt$Age_Group == min(split_propFilt$Age_Group)]
    split_propArray[split_propArray$age < min(split_propFilt$Age_Group),]$Prop_NSAS <- 
      split_propFilt$Prop_NSAS[split_propFilt$Age_Group == min(split_propFilt$Age_Group)]
    # all ages above maximum is set to maximum
    split_propArray[split_propArray$age > max(split_propFilt$Age_Group),]$Prop_WBSS <- 
      split_propFilt$Prop_WBSS[split_propFilt$Age_Group == max(split_propFilt$Age_Group)]
    split_propArray[split_propArray$age > max(split_propFilt$Age_Group),]$Prop_NSAS <- 
      split_propFilt$Prop_NSAS[split_propFilt$Age_Group == max(split_propFilt$Age_Group)]
    
    # loop on all the index ages (1-9+) with 9 as plus group
    for(idxAges in 1:nAges){
      
      # select age to be computed, make sure we combine ages for the plus group
      if(idxAges == nAges){
        # select ages as plut group
        ageSel <- unique(endTabNO[endTabNO$strata == strataCurrent,]$age)[unique(endTabNO[endTabNO$strata == strataCurrent,]$age) >= ageVec[idxAges]]
      }else{
        # select current age
        ageSel <- ageVec[idxAges]
      }
      
      indexTemp.NSAS <- array(0,dim=c(length(ageSel),1))
      indexTemp.WBSS <- array(0,dim=c(length(ageSel),1))
      # loop on all ages available in the StoX object
      for(idxAgeSel in 1:length(ageSel)){
        # filter superInd table
        idxTabNOFilt <- which(endTabNO$strata == strataCurrent & endTabNO$age %in% ageSel[idxAgeSel])
        
        if(ageSel[idxAgeSel] > max(split_propArray$age)){
          splitCurrent <- split_propArray[split_propArray$age == max(split_propArray$age),]
        }else{
          splitCurrent <- split_propArray[split_propArray$age == ageSel[idxAgeSel],]
        }
        
        endTabNO.NSAS$Abundance[idxTabNOFilt] <- endTabNO.NSAS$Abundance[idxTabNOFilt]*splitCurrent$Prop_NSAS
        endTabNO.WBSS$Abundance[idxTabNOFilt] <- endTabNO.WBSS$Abundance[idxTabNOFilt]*splitCurrent$Prop_WBSS
        endTabNO.NSAS$biomass_g[idxTabNOFilt] <- endTabNO.NSAS$biomass_g[idxTabNOFilt]*splitCurrent$Prop_NSAS
        endTabNO.WBSS$biomass_g[idxTabNOFilt] <- endTabNO.WBSS$biomass_g[idxTabNOFilt]*splitCurrent$Prop_WBSS
        
        # calculate split numbers
        indexTemp.NSAS[idxAgeSel] <- sum(endTabNO.NSAS$Abundance[idxTabNOFilt])
        indexTemp.WBSS[idxAgeSel] <- sum(endTabNO.WBSS$Abundance[idxTabNOFilt])
      }
      indexHERASComponent_NO.NSAS[idxStrata,idxAges] <- sum(indexTemp.NSAS)
      indexHERASComponent_NO.WBSS[idxStrata,idxAges] <- sum(indexTemp.WBSS)
    }
  }
  
  # assign stock ID
  endTabNO.NSAS$stock <- 'her-47d3'
  endTabNO.WBSS$stock <- 'her-3a22'
  
  # create superIndividual split table
  endTabNO <- rbind(endTabNO.NSAS,endTabNO.WBSS)
  
  #################################################################
  ## EU project
  #################################################################
  endTabEU <- fill_missing_species_v3(endTabEU,mode='Abundance')
  
  uniqueStrataEU <- unique(endTabEU$strata)

  # initialize array nStrata x nAges
  indexHERASComponent_EU.NSAS <- array(0, dim=c(length(uniqueStrataEU),nAges))
  indexHERASComponent_EU.WBSS <- array(0, dim=c(length(uniqueStrataEU),nAges))
  
  for(idxStrata in 1:length(uniqueStrataEU)){
    strataCurrent <- uniqueStrataEU[idxStrata]
    
    endTabEUFilt <- endTabEU[endTabEU$strata == strataCurrent,]
    
    # loop on all the index ages (1-9+) with 9 as plus group
    for(idxAges in 1:nAges){

      # select age to be computed, make sure we combine ages for the plus group
      if(idxAges == nAges){
        # select ages as plut grounp
        ageSel <- unique(endTabEUFilt$age)[unique(endTabEUFilt$age) >= ageVec[idxAges]]
      }else{
        # select current age
        ageSel <- ageVec[idxAges]
      }
      
      endTabFilt.NSAS  <- endTabEUFilt[ endTabEUFilt$age %in% ageSel &
                                        endTabEUFilt$stock == 'her-47d3',]
      endTabFilt.WBSS  <- endTabEUFilt[ endTabEUFilt$age %in% ageSel &
                                        endTabEUFilt$stock == 'her-3a22',]
      
      indexHERASComponent_EU.NSAS[idxStrata,idxAges] <- sum(endTabFilt.NSAS$Abundance)
      indexHERASComponent_EU.WBSS[idxStrata,idxAges] <- sum(endTabFilt.WBSS$Abundance)
    }
  }
  
  endTabAll <- rbind(endTabEU,endTabNO)
  #endTabAll <- endTabNO
  
  #################################################################
  ## compute proportion mature and weight at age
  #################################################################
  
  # agregate table
  endTabAgg <- aggregate(cbind(endTabAll$Abundance,endTabAll$biomass_g),
                          list(   endTabAll$stock,
                                  endTabAll$age,
                                  endTabAll$maturity),
                          sum)
  
  colnames(endTabAgg) <- c('stock','age','maturity','Abundance','biomass')
  endTabAgg$stock[endTabAgg$stock == 'her-47d3'] <- 'NSAS'
  endTabAgg$stock[endTabAgg$stock == 'her-3a22'] <- 'WBSS'
  
  #endTabAgg <- aggregate(cbind(endTabAll$Abundance,endTabAll$biomass_g),
  #                       list(   endTabAll$stock,
  #                               endTabAll$age,
  #                               endTabAll$strata,
  #                               endTabAll$maturity),
  #                       sum)
  #colnames(endTabAgg) <- c('stock','age','strata','maturity','Abundance','biomass')
  #endTabAgg$stock[endTabAgg$stock == 'her-47d3'] <- 'NSAS'
  #endTabAgg$stock[endTabAgg$stock == 'her-3a22'] <- 'WBSS'
  
  #endTabAgg[endTabAgg$strata == 141& endTabAgg$stock == 'WBSS'&endTabAgg$maturity == 'MAT',]
  
  # initialize weight at age vectors
  waMat <- as.data.frame(array(0, dim=c(nAges,3)))
  colnames(waMat) <- c('age','NSAS','WBSS')
  waMat$age <- ageVec
  
  # initialize proportion mature vectors
  propMat <- as.data.frame(array(0, dim=c(nAges,3)))
  colnames(propMat) <- c('age','NSAS','WBSS')
  propMat$age <- ageVec
  
  for(idxAge in ageVec){
    for(idxStock in c('NSAS','WBSS')){
      if(idxAges == max(ageVec)){
        # proportion mature
        propMat[propMat$age == idxAge,idxStock] <- 
          sum(endTabAgg[endTabAgg$age >= idxAge & endTabAgg$stock == idxStock & endTabAgg$maturity == 'MAT',]$Abundance)/
          sum(endTabAgg[endTabAgg$age >= idxAge & endTabAgg$stock == idxStock,]$Abundance)
        
        # weight at age
        waMat[waMat$age == idxAge,idxStock] <- 
          sum(endTabAgg[endTabAgg$age >= idxAge & endTabAgg$stock == idxStock,]$biomass)/
          sum(endTabAgg[endTabAgg$age >= idxAge & endTabAgg$stock == idxStock,]$Abundance)
        
      }else{
        # proportion mature
        propMat[propMat$age == idxAge,idxStock] <- 
          sum(endTabAgg[endTabAgg$age == idxAge & endTabAgg$stock == idxStock & endTabAgg$maturity == 'MAT',]$Abundance)/
          sum(endTabAgg[endTabAgg$age == idxAge & endTabAgg$stock == idxStock,]$Abundance)
        
        # weight at age
        waMat[waMat$age == idxAge,idxStock] <- 
          sum(endTabAgg[endTabAgg$age == idxAge & endTabAgg$stock == idxStock,]$biomass)/
          sum(endTabAgg[endTabAgg$age == idxAge & endTabAgg$stock == idxStock,]$Abundance)
      }
    }
  }

  #propMat[is.na(propMat)] <- 0
  
  #################################################################
  ## tidying up before returning outputs
  #################################################################
  
  indexHERASComponent_EU.NSAS <- as.data.frame(indexHERASComponent_EU.NSAS,
                                               row.names=uniqueStrataEU)
  colnames(indexHERASComponent_EU.NSAS) <- ageVec
  
  indexHERASComponent_EU.WBSS <- as.data.frame(indexHERASComponent_EU.WBSS,
                                               row.names=uniqueStrataEU)
  colnames(indexHERASComponent_EU.WBSS) <- ageVec
  
  indexHERASComponent_NO.NSAS <- as.data.frame(indexHERASComponent_NO.NSAS,
                                               row.names=uniqueStrataNO)
  colnames(indexHERASComponent_NO.NSAS) <- ageVec
  
  indexHERASComponent_NO.WBSS <- as.data.frame(indexHERASComponent_NO.WBSS,
                                               row.names=uniqueStrataNO)
  colnames(indexHERASComponent_NO.WBSS) <- ageVec
  
  return(list(NSAS_EU = indexHERASComponent_EU.NSAS,
              WBSS_EU = indexHERASComponent_EU.WBSS,
              NSAS_NO = indexHERASComponent_NO.NSAS,
              WBSS_NO = indexHERASComponent_NO.WBSS,
              propM   = propMat,
              wa      = waMat))
  
}