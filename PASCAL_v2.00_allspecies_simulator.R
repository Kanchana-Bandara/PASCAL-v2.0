Simulator_1YGT_3GV1 <- function(){
      ################################################################################
      #Life-cycle simulator function :: 1-year life cycle :: 3RD GENERATION :: V.1.00#
      ################################################################################
      
      #Initial variables
      #+++++++++++++++++
      MLD <- MixedLayerDepth
      S <- SpeciesID
      I <- CarrierID
      phi <- Latitude
      k <- 8.62e-5
      T0 <- 273
      E_D <- .subset2(SSGDParameters, S)[13]
      a <- 0.6
      FCC <- ChlC
      FCMax <- MaxFA * FCC
      
      MC1 <- AMC_L
      MC2 <- AMC_H
      MC3 <- DMC
      p1 <- 40.29524                                                                #Exponent for starvation risk model
      p2 <- 0.20                                                                    #Coefficient of allometric scaling of p6
      p3 <- -0.175                                                                  #Exponent of allometric scaling of p6
      p4 <- 1                                                                       #Horizontal asymptote of food-limitation model
      p5 <- 0                                                                       #Value of p4 when food concentration = 0
      p6 <- NULL                                                                    #Numerically estimated as p2 * CurrentCWeight^p3
      WCF <- (1 / (40 / 100))
      
      StartPing <- STime[I]
      ISP <- isp[I]
      GAP <- gap[I]
      SDP <- sdp[I]
      SAP <- sap[I]
      LCLP <- lclp[I]
      StartDepth <- 1
      StartCWeight <- .subset2(SSGDParameters, S)[1]
      CurrentPing <- StartPing
      CurrentDepth <- StartDepth
      CurrentCWeight <- StartCWeight
      K <- VPRScalar
      CurrentNVPRisk <- NVPRScalar
      Suvivorship <- 1
      LifeState <- as.factor("ALIVE")
      CessasionStage <- NULL
      J <- 0
      CatCWeight <- 0
      STTreshold <- StarvationTolerance
      CriticalMoltingMasses <- CM_Lower + ((CM_Upper - CM_Lower) * LCLP)
      CMM <- CriticalMoltingMasses
      CurrentSWeight <- 0
      CurrentRWeight <- 0
      MaxCWeight <- NULL
      Fecundity <- 0
      EggWeight <- .subset2(SSGDParameters, S)[1]
      CB <- 0
      CBE <- 0
      CurrentFitness <- 0
      CessasionHorizon_T <- LCBreak
      EPState <- 0
      EPStartPing <- NULL
      
      IC_m <- SSGDParameters[2, S]
      IE_m <- SSGDParameters[3, S]
      IC_t <- SSGDParameters[4, S]
      IE_t <- SSGDParameters[5, S]
      MC_m <- SSGDParameters[6, S]
      ME_m <- SSGDParameters[7, S]
      MC_t <- SSGDParameters[8, S]
      ME_t <- SSGDParameters[9, S]
      dMC_m <- 6.8177e-4
      dME_m <- 0.7503
      dMC_t <- 1.1216
      dME_t <- 0.0875
      
      IterationID <- 0
      DTTrack <- NULL
      
      #DiapauseLog <- matrix(rep(NA, 10 * 2), ncol = 10)
      #colnames(DiapauseLog) <- c("DEnt", "DStg", "DIni", "DIniCW", "DIniSW", "ODepth", "DTer", "DTerCW", "DTerSW", "DStatus")
      #rownames(DiapauseLog) <- c("Year-1", "Year-2")
      #DepthLog <- rep(NA, 8760 * 2)
      #WeightLog <- rep(NA, 8760 * 2)
      #ReserveLog <- rep(NA, 8760 * 2)
      #SSLog <- rep(NA, 8760 * 2)
      #FecundityLog <- rep(NA, 8760 * 2)
      #DTLog <- rep(NA, 12)
      #ETLog <- rep(NA, 12)
      
      #Embroyonic development
      D0 <- .subset2(SSGDParameters, S)[10]
      
      repeat{
            #Timecoding and basic estimations
            IterationID <- IterationID + 1
            if(IterationID == 1){ CurrentPing <- CurrentPing }else{ CurrentPing <- CurrentPing + 1 }
            if(IterationID == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- ceiling(MLD[CurrentPing] * rbeta(1, 1, 5)) }
            CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
            if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
            
            #Physical submodel :: extracting time specific vectors
            CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
            CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
            
            #Development time estimation
            CurrentTemp <- CurrentTempRange[CurrentDepth] + 273.15
            CurrentDT <- round(D0 * exp((-E_D * (CurrentTemp - T0)) / (k * CurrentTemp * T0)) * 24, 2)
            DTTrack <- append(DTTrack, CurrentDT)
            
            #Degrowth estimation
            CurrentTemp <- CurrentTempRange[CurrentDepth]
            MRate_m <- MC_m * CurrentCWeight^ME_m
            MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
            CurrentMRate <- MRate_m * MRate_t_scalar * 0.5
            CurrentCWeight <- CurrentCWeight - CurrentMRate
            
            #Mortality risk estimation
            CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
            KAdj <- K_base[ceiling(CurrentCWeight)] * K
            CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
            if(CurrentCWeight <= 0){ CurrentSRisk <- 1 }else{ CurrentSRisk <- 0 }
            Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
            
            #Tracking and high-resolution logging (inactive during simulation)
            #-----------------------------------------------------------------
            #DepthLog[CurrentPing] <- CurrentDepth 
            #WeightLog[CurrentPing] <- CurrentCWeight
            #SSLog[CurrentPing] <- Suvivorship
            #-----------------------------------------------------------------
            
            #Death or development time based breaking condition
            if(Suvivorship < 1e-2){
                  LifeState <- as.factor("DEAD")
                  CessasionStage <- J
                  IterationID <- 0
                  DTTrack <- NULL
                  break
            }else if(IterationID >= gm_mean(DTTrack)){
                  J <- J + 1
                  #DTLog[J] <- IterationID
                  #ETLog[J] <- CurrentPing
                  IterationID <- 0
                  DTTrack <- NULL
                  break
            }else{
                  #proceed iteration
            }
      }
      
      #Non-feeding nauplii (NI-NII)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            D0SSV <- c(.subset2(SSGDParameters, S)[11], .subset2(SSGDParameters, S)[12])
            
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  D0 <- D0SSV[J]
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
                  CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
                  CurrentIRadRange <- .subset2(Irradiance, CPRemap)[1:MaxDepth]
                  
                  #Distance and range estimates
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
                  
                  #Diel depth selection
                  if(any(CurrentIRadRange[SR] < ISPAdj)){
                        CurrentDepth <- SR[which(CurrentIRadRange[SR] < ISPAdj)]                                                                                        
                        CurrentTemp <- CurrentTempRange[CurrentDepth] + 273.15                                                                                                   
                        CurrentDT <- round(D0 * exp((-E_D * (CurrentTemp - T0)) / (k * CurrentTemp * T0)) * 24, 2)                                                      
                        CurrentDepth <- CurrentDepth[which(CurrentDT == min(CurrentDT))]
                        if(length(CurrentDepth) == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- min(CurrentDepth) }
                  }else{
                        CurrentDepth <- SR[length(SR)]                                                                                                                  
                  }
                  
                  #Estimations based on the current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Development time estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth] + 273.15
                  CurrentDT <- round(D0 * exp((-E_D * (CurrentTemp - T0)) / (k * CurrentTemp * T0)) * 24, 2)
                  DTTrack <- append(DTTrack, CurrentDT)
                  
                  #Degrowth estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth]
                  MRate_m <- MC_m * CurrentCWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar * 0.5
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentCWeight <- CurrentCWeight - CurrentMRate
                  
                  #Mortality risk estimation
                  CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  if(CurrentCWeight <= 0){ CurrentSRisk <- 1 }else{ CurrentSRisk <- 0 }
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #-----------------------------------------------------------------
                  
                  #Death or development time based breaking condition
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL                                                                                                                                       
                        break
                  }else if(IterationID >= gm_mean(DTTrack)){
                        J <- J + 1
                        #DTLog[J] <- IterationID
                        #ETLog[J] <- CurrentPing
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL
                        if(J >= 3){
                              break
                        }else{
                              #Do not break
                        }
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      #Feeding stages with no energetic reserve (NIII - CIII)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            MaxCWeight <- CurrentCWeight
            
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
                  CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
                  CurrentIRadRange <- .subset2(Irradiance, CPRemap)[1:MaxDepth]
                  CurrentFARange <- .subset2(FA, CPRemap)[1:MaxDepth]
                  
                  #Distance and range estimates
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
                  
                  #Diel-depth selection
                  if(any(CurrentIRadRange[SR] < ISPAdj)){
                        CurrentDepth <- SR[which(CurrentIRadRange[SR] < ISPAdj)]                                                                                        
                        CurrentTemp <- CurrentTempRange[CurrentDepth]
                        CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                        p6 <- p2 * CurrentCWeight^p3
                        CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                        IRate_m <- IC_m * CurrentCWeight^IE_m
                        IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                        CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                        CurrentGPotential <- a * CurrentIRate
                        CurrentDepth <- CurrentDepth[which(CurrentGPotential == max(CurrentGPotential))]                                                                                
                        if(length(CurrentDepth) == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- min(CurrentDepth) }
                  }else{
                        CurrentDepth <- SR[length(SR)]
                  }
                  
                  #Estimations based on current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Growth estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth]
                  CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                  p6 <- p2 * CurrentCWeight^p3
                  CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                  IRate_m <- IC_m * CurrentCWeight^IE_m
                  IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                  CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                  MRate_m <- MC_m * CurrentCWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentGPotential <- (a * CurrentIRate) - CurrentMRate
                  CurrentCWeight <- CurrentCWeight + CurrentGPotential
                  
                  if(CurrentCWeight > MaxCWeight){ MaxCWeight <- CurrentCWeight }else{ MaxCWeight <- MaxCWeight }
                  
                  #Mortality risk estimation :: visual predation, non-visual predation and starvation risks (death at 50% structural mass loss :: Chossat, 1843)
                  CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  CatCWeight <- (MaxCWeight - CurrentCWeight) / MaxCWeight
                  
                  if(CatCWeight <= STTreshold){
                        CurrentSRisk <- 0
                  }else if(CatCWeight > STTreshold & CatCWeight <= 0.50){
                        CurrentSRisk <- 1e-7 * exp(p1 * CatCWeight)
                  }else{
                        CurrentSRisk <- 1
                  }
                  
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #-----------------------------------------------------------------
                  
                  #Death or critical molting mass based breaking condition
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL                                                                                                                                       
                        break
                  }else if(CurrentCWeight >= CMM[J + 1]){
                        J <- J + 1
                        #DTLog[J] <- IterationID
                        #ETLog[J] <- CurrentPing
                        IterationID <- 0
                        if(J >= 10){
                              break
                        }else{
                              #Don't break
                        }
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      #Feeding stages with an energetic reserve (CIV-CV)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
                  CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
                  CurrentIRadRange <- .subset2(Irradiance, CPRemap)[1:MaxDepth]
                  CurrentFARange <- .subset2(FA, CPRemap)[1:MaxDepth]
                  
                  #Distance and range estimates
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
                  
                  #Diel-depth selection
                  if(any(CurrentIRadRange[SR] < ISPAdj)){
                        CurrentDepth <- SR[which(CurrentIRadRange[SR] < ISPAdj)]                                                                                        
                        CurrentTemp <- CurrentTempRange[CurrentDepth]
                        CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                        p6 <- p2 * CurrentCWeight^p3
                        CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                        IRate_m <- IC_m * CurrentCWeight^IE_m
                        IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                        CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                        CurrentGPotential <- a * CurrentIRate
                        CurrentDepth <- CurrentDepth[which(CurrentGPotential == max(CurrentGPotential))]                                                                                
                        if(length(CurrentDepth) == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- min(CurrentDepth) }
                  }else{
                        CurrentDepth <- SR[length(SR)]
                  }
                  
                  #Estimations based on current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Growth estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth]
                  CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                  p6 <- p2 * CurrentCWeight^p3
                  CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                  IRate_m <- IC_m * CurrentCWeight^IE_m
                  IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                  CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                  CurrentTWeight <- CurrentCWeight + CurrentSWeight
                  MRate_m <- MC_m * CurrentTWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentGPotential <- (a * CurrentIRate) - CurrentMRate
                  
                  if(CurrentGPotential >= 0){
                        if(CurrentCWeight >= max(CMM)){
                              CurrentCWeight <- CurrentCWeight
                              CurrentSWeight <- CurrentSWeight + (CurrentGPotential * SVol[ceiling(CurrentCWeight)])
                        }else{
                              MaxSAllocation <- CurrentGPotential * GAP * SVol[ceiling(CurrentCWeight)]
                              CurrentSWeight <- CurrentSWeight + MaxSAllocation
                              CurrentCWeight <- CurrentCWeight + (CurrentGPotential - MaxSAllocation)
                        }
                  }else{
                        if(CurrentSWeight >= abs(CurrentGPotential)){
                              CurrentCWeight <- CurrentCWeight
                              CurrentSWeight <- CurrentSWeight - abs(CurrentGPotential)
                        }else{
                              CurrentCWeight <- CurrentCWeight + (CurrentGPotential + CurrentSWeight)
                              CurrentSWeight <- 0
                        }
                  }
                  
                  if(CurrentCWeight > MaxCWeight){ MaxCWeight <- CurrentCWeight }else{ MaxCWeight <- MaxCWeight }
                  
                  #Mortality risk estimation :: visual predation, non-visual predation and starvation risks (death at 50% structural mass loss :: Chossat, 1843)
                  CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  CatCWeight <- (MaxCWeight - CurrentCWeight) / MaxCWeight
                  
                  if(CatCWeight <= STTreshold){
                        CurrentSRisk <- 0
                  }else if(CatCWeight > STTreshold & CatCWeight <= 0.50){
                        CurrentSRisk <- 1e-7 * exp(p1 * CatCWeight)
                  }else{
                        CurrentSRisk <- 1
                  }
                  
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #ReserveLog[CurrentPing] <- CurrentSWeight
                  #-----------------------------------------------------------------
                  
                  #Death or critical molting mass based breaking condition or diapause entry
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL                                                                                                                                       
                        break
                  }else if(J == 10){
                        if(CurrentCWeight >= CMM[J + 1]){
                              J <- J + 1
                              #DTLog[J] <- IterationID
                              #ETLog[J] <- CurrentPing
                              IterationID <- 0
                        }else if((CurrentSWeight / CurrentCWeight) >= SDP){
                              #DiapauseLog[1, 1] <- CurrentPing
                              #DiapauseLog[1, 2] <- J
                              break
                        }else{
                              #Proceed iteration
                        }
                  }else if((CurrentSWeight / CurrentCWeight) >= SDP){
                        #DiapauseLog[1, 1] <- CurrentPing
                        #DiapauseLog[1, 2] <- J
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      #Diapause entry (CIV, CV)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            #Overwintering depth selection
            OD_T <- max(MixedLayerDepth)
            PotentialOD_I <- .subset2(Irradiance, CPRemap)[MinDepth:MaxDepth]
            ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
            OD_I <- min(which(PotentialOD_I <= ISPAdj))
            if(OD_I >= OD_T){ OD_Min <- OD_I }else{ OD_Min <- OD_T }
            
            repeat{
                  OD <- round(rnorm(1, mean = mean(OD_Min:MaxDepth), sd = (OD_Min / 10)), 0)
                  if(OD <= MaxDepth){
                        break
                  }
            }
            
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  #Search for overwintering depth
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  
                  if(length(which(SR == OD)) == 1){
                        CurrentDepth <- OD
                  }else{
                        if(max(SR) > OD){ CurrentDepth <- SR[1] }else{ CurrentDepth <- SR[length(SR)] }
                  }
                  
                  #Estimations based on current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTemp <- .subset2(Temperature, CPRemap)[CurrentDepth]
                  CurrentIVPR <- .subset2(PredationRisk, CPRemap)[CurrentDepth]
                  
                  #Degrowth estimation
                  CurrentTWeight <- CurrentCWeight + CurrentSWeight
                  MRate_m <- MC_m * CurrentTWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentGPotential <- 0 - CurrentMRate
                  
                  if(CurrentSWeight >= abs(CurrentGPotential)){
                        CurrentCWeight <- CurrentCWeight
                        CurrentSWeight <- CurrentSWeight - abs(CurrentGPotential)
                  }else{
                        CurrentCWeight <- CurrentCWeight + (CurrentGPotential + CurrentSWeight)
                        CurrentSWeight <- 0
                  }
                  
                  if(CurrentCWeight > MaxCWeight){ MaxCWeight <- CurrentCWeight }else{ MaxCWeight <- MaxCWeight }
                  
                  #Mortality risk estimation :: visual predation, non-visual predation and starvation risks (death at 50% structural mass loss :: Chossat, 1843)
                  CurrentVPRisk <- CurrentIVPR
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  CatCWeight <- (MaxCWeight - CurrentCWeight) / MaxCWeight
                  
                  if(CatCWeight <= STTreshold){
                        CurrentSRisk <- 0
                  }else if(CatCWeight > STTreshold & CatCWeight <= 0.50){
                        CurrentSRisk <- 1e-7 * exp(p1 * CatCWeight)
                  }else{
                        CurrentSRisk <- 1
                  }
                  
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #ReserveLog[CurrentPing] <- CurrentSWeight
                  #-----------------------------------------------------------------
                  
                  #Death or diapause progression
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        break
                  }else if(CurrentDepth == OD){
                        DEntry[I] <<- CurrentPing
                        DStage[I] <<- J
                        DCWeight[I] <<- CurrentCWeight
                        DEntrySWeight[I] <<- CurrentSWeight
                        DDepth[I] <<- CurrentDepth
                        #DiapauseLog[1, 3] <- CurrentPing
                        #DiapauseLog[1, 4] <- CurrentCWeight
                        #DiapauseLog[1, 5] <- CurrentSWeight
                        #DiapauseLog[1, 6] <- CurrentDepth
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      if(DiapauseMachineCode == "OFF"){
            #Diapause progression (CIV, CV) :: Full-length code
            if(LifeState == "DEAD"){
                  #Do nothing and proceed
            }else{
                  TempRange <- Temperature[OD, ]
                  IVPRRange <- PredationRisk[OD, ]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  InitialSWeight <- CurrentSWeight
                  
                  repeat{
                        IterationID <- IterationID + 1
                        CurrentPing <- CurrentPing + 1
                        CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                        if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                        
                        CurrentTemp <- TempRange[CPRemap]
                        CurrentTWeight <- CurrentCWeight + CurrentSWeight
                        MRate_m <- dMC_m * CurrentTWeight^dME_m
                        MRate_t_scalar <- dMC_t * exp(CurrentTemp * dME_t)
                        CurrentMRate <- MRate_m * MRate_t_scalar * MC3
                        CurrentSWeight <- CurrentSWeight - CurrentMRate
                        CurrentCWeight <- CurrentCWeight
                        
                        CurrentVPRisk <- IVPRRange[CPRemap]
                        CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                        Suvivorship <- Suvivorship * (1 - CurrentPRisk)
                        
                        #Tracking and high-resolution logging (inactive during simulation)
                        #-----------------------------------------------------------------
                        #DepthLog[CurrentPing] <- CurrentDepth 
                        #WeightLog[CurrentPing] <- CurrentCWeight
                        #SSLog[CurrentPing] <- Suvivorship
                        #ReserveLog[CurrentPing] <- CurrentSWeight
                        #-----------------------------------------------------------------
                        
                        if(Suvivorship < 1e-2){
                              #DiapauseLog[1, 7] <- CurrentPing
                              #DiapauseLog[1, 8] <- CurrentCWeight
                              #DiapauseLog[1, 9] <- CurrentSWeight
                              #DiapauseLog[1, 10] <- 0
                              break
                        }else if(CurrentSWeight <= InitialSWeight * (1 - SAP)){
                              DExit[I] <<- CurrentPing
                              DExitSWeight[I] <<- CurrentSWeight
                              #DiapauseLog[1, 7] <- CurrentPing
                              #DiapauseLog[1, 8] <- CurrentCWeight
                              #DiapauseLog[1, 9] <- CurrentSWeight
                              #DiapauseLog[1, 10] <- 1
                              break
                        }else{
                              #Proceed iteration
                        }
                  }
            }
      }else{
            #Diapause progression (CIV, CV) :: Machine code :: Logging not possible
            if(LifeState == "DEAD"){
                  #Do nothing and proceed
            }else{
                  CurrentPing <- CurrentPing + 1
                  IterationID <- IterationID + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  CurrentTemp <- .subset2(Temperature, CPRemap)[CurrentDepth]
                  CurrentIVPR <- .subset2(PredationRisk, CPRemap)[CurrentDepth]
                  
                  MRate_t_scalar <- dMC_t * exp(CurrentTemp * dME_t)
                  CurrentTWeight <- CurrentCWeight + CurrentSWeight
                  MRate_m_ceiling <- dMC_m * CurrentTWeight^dME_m
                  MRate_ceiling <- MRate_m_ceiling * MRate_t_scalar * MC3
                  ProjectedSWeight <- CurrentSWeight - (CurrentSWeight * SAP)
                  ProjectedTWeight <- CurrentCWeight + ProjectedSWeight
                  MRate_m_floor <- dMC_m * ProjectedTWeight^dME_m
                  Mrate_floor <- MRate_m_floor * MRate_t_scalar * MC3
                  PredictedMRate <- (MRate_ceiling + Mrate_floor) / 2
                  DiapauseDuration <- round((CurrentSWeight * SAP) / PredictedMRate, 0)
                  CurrentSWeight <- ProjectedSWeight
                  CurrentCWeight <- CurrentCWeight
                  
                  CurrentVPRisk <- CurrentIVPR
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  DiapauseSuvivorship <- (1 - CurrentPRisk)^DiapauseDuration
                  Suvivorship <- Suvivorship * DiapauseSuvivorship
                  
                  CurrentPing <- CurrentPing + DiapauseDuration
                  IterationID <- IterationID + DiapauseDuration
                  
                  if(Suvivorship < 1e-2){
                        #DiapauseLog[1, 7] <- CurrentPing
                        #DiapauseLog[1, 8] <- CurrentCWeight
                        #DiapauseLog[1, 9] <- CurrentSWeight
                        #DiapauseLog[1, 10] <- 0
                  }else{
                        DExit[I] <<- CurrentPing
                        DExitSWeight[I] <<- CurrentSWeight
                        #DiapauseLog[1, 7] <- CurrentPing
                        #DiapauseLog[1, 8] <- CurrentCWeight
                        #DiapauseLog[1, 9] <- CurrentSWeight
                        #DiapauseLog[1, 10] <- 1
                  }
            }
      }
      
      #Diapause exit (CIV, CV)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
                  CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
                  CurrentIRadRange <- .subset2(Irradiance, CPRemap)[1:MaxDepth]
                  CurrentFARange <- .subset2(FA, CPRemap)[1:MaxDepth]
                  
                  #Distance and range estimates
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
                  
                  #Diel-depth selection
                  if(any(CurrentIRadRange[SR] < ISPAdj)){
                        CurrentDepth <- SR[which(CurrentIRadRange[SR] < ISPAdj)]                                                                                        
                        CurrentTemp <- CurrentTempRange[CurrentDepth]
                        CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                        p6 <- p2 * CurrentCWeight^p3
                        CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                        IRate_m <- IC_m * CurrentCWeight^IE_m
                        IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                        CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                        CurrentGPotential <- a * CurrentIRate
                        CurrentDepth <- CurrentDepth[which(CurrentGPotential == max(CurrentGPotential))]                                                                                
                        if(length(CurrentDepth) == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- min(CurrentDepth) }
                  }else{
                        CurrentDepth <- SR[length(SR)]
                  }
                  
                  #Estimations based on current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Growth estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth]
                  CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                  p6 <- p2 * CurrentCWeight^p3
                  CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                  IRate_m <- IC_m * CurrentCWeight^IE_m
                  IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                  CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                  CurrentTWeight <- CurrentCWeight + CurrentSWeight
                  MRate_m <- MC_m * CurrentTWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentGPotential <- (a * CurrentIRate) - CurrentMRate
                  
                  if(CurrentGPotential >= 0){
                        if(CurrentCWeight >= max(CMM)){
                              CurrentCWeight <- CurrentCWeight
                              CurrentSWeight <- CurrentSWeight
                        }else{
                              CurrentCWeight <- CurrentCWeight + CurrentGPotential
                              CurrentSWeight <- CurrentSWeight
                        }
                  }else{
                        if(CurrentSWeight >= abs(CurrentGPotential)){
                              CurrentCWeight <- CurrentCWeight
                              CurrentSWeight <- CurrentSWeight - abs(CurrentGPotential)
                        }else{
                              CurrentCWeight <- CurrentCWeight + (CurrentGPotential + CurrentSWeight)
                              CurrentSWeight <- 0
                        }
                  }
                  
                  if(CurrentCWeight > MaxCWeight){ MaxCWeight <- CurrentCWeight }else{ MaxCWeight <- MaxCWeight }
                  
                  #Mortality risk estimation :: visual predation, non-visual predation and starvation risks (death at 50% structural mass loss :: Chossat, 1843)
                  CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  CatCWeight <- (MaxCWeight - CurrentCWeight) / MaxCWeight
                  
                  if(CatCWeight <= STTreshold){
                        CurrentSRisk <- 0
                  }else if(CatCWeight > STTreshold & CatCWeight <= 0.50){
                        CurrentSRisk <- 1e-7 * exp(p1 * CatCWeight)
                  }else{
                        CurrentSRisk <- 1
                  }
                  
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #ReserveLog[CurrentPing] <- CurrentSWeight
                  #-----------------------------------------------------------------
                  
                  #Death or critical molting mass based breaking condition
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        break
                  }else if(CurrentCWeight >= CMM[J + 1]){
                        J <- J + 1
                        #DTLog[J] <- IterationID
                        #ETLog[J] <- CurrentPing
                        IterationID <- 0
                        if(J == 12){
                              break
                        }else{
                              #Dont break
                        }
                  }else{
                        #Proceed iteration
                  }
            }
      }
      
      #Reproductive stage (CVI-F)
      if(LifeState == "DEAD"){
            #Do nothing and proceed
      }else{
            repeat{
                  #Timecoding and basic estimations
                  IterationID <- IterationID + 1
                  CurrentPing <- CurrentPing + 1
                  CPRemap <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                  if(CPRemap <= 0){ CPRemap <- 1 }else{ CPRemap <- CPRemap }
                  
                  #Physical submodel :: extracting time specific vectors
                  CurrentTempRange <- .subset2(Temperature, CPRemap)[1:MaxDepth]
                  CurrentIVPRRange <- .subset2(PredationRisk, CPRemap)[1:MaxDepth]
                  CurrentIRadRange <- .subset2(Irradiance, CPRemap)[1:MaxDepth]
                  CurrentFARange <- .subset2(FA, CPRemap)[1:MaxDepth]
                  
                  #Distance and range estimates
                  MaxDistance <- round(8.0116 * (CurrentCWeight^0.4531), 0)
                  if((CurrentDepth + MaxDistance) > MaxDepth){ MFloor <- MaxDepth }else{ MFloor <- CurrentDepth + MaxDistance }
                  if((CurrentDepth - MaxDistance) < MinDepth){ MCeiling <- MinDepth }else{ MCeiling <- CurrentDepth - MaxDistance }
                  SR <- MCeiling:MFloor
                  InitialDepth <- CurrentDepth
                  ISPAdj <- ISP_base_rc[ceiling(CurrentCWeight)] * ISP
                  
                  #Diel-depth selection
                  if(any(CurrentIRadRange[SR] < ISPAdj)){
                        CurrentDepth <- SR[which(CurrentIRadRange[SR] < ISPAdj)]                                                                                        
                        CurrentTemp <- CurrentTempRange[CurrentDepth]
                        CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                        p6 <- p2 * CurrentCWeight^p3
                        CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                        IRate_m <- IC_m * CurrentCWeight^IE_m
                        IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                        CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                        CurrentGPotential <- a * CurrentIRate
                        CurrentDepth <- CurrentDepth[which(CurrentGPotential == max(CurrentGPotential))]                                                                                
                        if(length(CurrentDepth) == 1){ CurrentDepth <- CurrentDepth }else{ CurrentDepth <- min(CurrentDepth) }
                  }else{
                        CurrentDepth <- SR[length(SR)]
                  }
                  
                  #Estimations based on current and initial depths
                  MovingTime <- abs(CurrentDepth - InitialDepth) / MaxDistance
                  IdleTime <- 1 - MovingTime
                  VMTrajectory <- InitialDepth:CurrentDepth
                  VMSML <- round(length(which(VMTrajectory <= MLD[CPRemap])) / length(VMTrajectory), 2)
                  VMHGZ <- 1 - VMSML
                  
                  #Growth estimation
                  CurrentTemp <- CurrentTempRange[CurrentDepth]
                  CurrentFA <- CurrentFARange[CurrentDepth] * FCC
                  p6 <- p2 * CurrentCWeight^p3
                  CurrentFL <- p4 + (p5 - p4) * exp(-p6 * CurrentFA)
                  IRate_m <- IC_m * CurrentCWeight^IE_m
                  IRate_t_scalar <-  IC_t * exp(CurrentTemp * IE_t)
                  CurrentIRate <-IRate_m * IRate_t_scalar * CurrentFL
                  CurrentTWeight <- CurrentCWeight + CurrentSWeight
                  MRate_m <- MC_m * CurrentTWeight^ME_m
                  MRate_t_scalar <- MC_t * exp(CurrentTemp * ME_t)
                  CurrentBMRate <- MRate_m * MRate_t_scalar
                  CurrentAMRate <- (CurrentBMRate * MC1 * MovingTime * VMHGZ) + (CurrentBMRate * MC2 * MovingTime * VMSML)
                  CurrentMRate <- CurrentBMRate + CurrentAMRate
                  CurrentGPotential <- (a * CurrentIRate) - CurrentMRate
                  
                  #Reproduction
                  MaxIRate_m <- IC_m * CurrentCWeight^IE_m
                  IRate_t_scalar <- IC_t * exp(CurrentTemp * IE_t)
                  MaxIRate <- MaxIRate_m * IRate_t_scalar
                  MaxGPotential <- a * MaxIRate
                  
                  if(CurrentGPotential >= 0){
                        if(S == 2){
                              CurrentCWeight <- CurrentCWeight
                              CurrentRWeight <- CurrentRWeight + CurrentGPotential
                              CurrentSWeight <- CurrentSWeight
                              CB <- 0
                        }else if(S == 3){
                              if(CurrentSWeight >= MaxGPotential){
                                    if(MaxGPotential > CurrentGPotential){
                                          CurrentCWeight <- CurrentCWeight
                                          MaxRMixAllocation <- MaxGPotential - CurrentGPotential
                                          CurrentRWeight <- CurrentRWeight + CurrentGPotential + MaxRMixAllocation
                                          CurrentSWeight <- CurrentSWeight - MaxRMixAllocation
                                          CB <- 1
                                    }else{
                                          CurrentCWeight <- CurrentCWeight
                                          CurrentRWeight <- CurrentRWeight + CurrentGPotential
                                          CurrentSWeight <- CurrentSWeight
                                          CB <- 0
                                    }
                              }else{
                                    CurrentCWeight <- CurrentCWeight
                                    CurrentRWeight <- CurrentRWeight + CurrentGPotential
                                    CurrentSWeight <- CurrentSWeight
                                    CB <- 0
                              }
                        }else{
                              if(CurrentSWeight >= MaxGPotential){
                                    CurrentCWeight <- CurrentCWeight
                                    CurrentSWeight <- CurrentSWeight - MaxGPotential
                                    CurrentRWeight <- CurrentRWeight + MaxGPotential
                                    CB <- 1
                              }else{
                                    CurrentCWeight <- CurrentCWeight
                                    CurrentRWeight <- CurrentRWeight + 0
                                    CurrentSWeight <- CurrentSWeight
                                    CB <- 0
                              }
                        }
                  }else{
                        if(CurrentSWeight >= MaxGPotential){
                              if(S == 2){
                                    CurrentCWeight <- CurrentCWeight
                                    CurrentSWeight <- CurrentSWeight - abs(CurrentGPotential)
                                    CurrentRWeight <- CurrentRWeight
                                    CB <- 0
                              }else{
                                    CurrentCWeight <- CurrentCWeight
                                    CurrentSWeight <- CurrentSWeight - MaxGPotential
                                    CurrentRWeight <- CurrentRWeight + (MaxGPotential - abs(CurrentGPotential))
                                    if(CurrentRWeight >= 0){ CB <- 1 }else{ CB <- 0 }
                              }
                        }else{
                              CurrentCWeight <- CurrentCWeight + (CurrentGPotential + CurrentSWeight)
                              CurrentSWeight <- 0
                              CurrentRWeight <- CurrentRWeight
                              CB <- 0
                        }
                  }
                  
                  if(CurrentRWeight >= EggWeight){
                        NEggs <- floor(CurrentRWeight / EggWeight)
                        Fecundity <- Fecundity + NEggs
                        CurrentRWeight <- CurrentRWeight %% EggWeight
                        if(CB == 1){ CBE <- CBE + NEggs }else{ CBE <- CBE }
                  }else{
                        NEggs <- 0
                        Fecundity <- Fecundity + NEggs
                        CurrentRWeight <- CurrentRWeight
                        CBE <- CBE
                  }
                  
                  if(Fecundity > 0 & EPState == 0){
                        EPState <- 1
                        EPStartPing <- CurrentPing
                  }
                  
                  #Mortality risk estimation :: visual predation, non-visual predation and starvation risks (death at 50% structural mass loss :: Chossat, 1843)
                  CurrentVPRisk <- CurrentIVPRRange[CurrentDepth]
                  KAdj <- K_base[ceiling(CurrentCWeight)] * K
                  CurrentPRisk <- (CurrentVPRisk + CurrentNVPRisk) * KAdj
                  CatCWeight <- (MaxCWeight - CurrentCWeight) / MaxCWeight
                  
                  if(CatCWeight <= STTreshold){
                        CurrentSRisk <- 0
                  }else if(CatCWeight > STTreshold & CatCWeight <= 0.50){
                        CurrentSRisk <- 1e-7 * exp(p1 * CatCWeight)
                  }else{
                        CurrentSRisk <- 1
                  }
                  
                  Suvivorship <- Suvivorship * (1 - (CurrentPRisk + CurrentSRisk))
                  CurrentFitness <- CurrentFitness + (NEggs * Suvivorship)
                  
                  #Tracking and high-resolution logging (inactive during simulation)
                  #-----------------------------------------------------------------
                  #DepthLog[CurrentPing] <- CurrentDepth 
                  #WeightLog[CurrentPing] <- CurrentCWeight
                  #SSLog[CurrentPing] <- Suvivorship
                  #ReserveLog[CurrentPing] <- CurrentSWeight
                  #FecundityLog[CurrentPing] <- Fecundity
                  #-----------------------------------------------------------------
                  
                  #Death or force stop and fitness weighting
                  if(Suvivorship < 1e-2){
                        LifeState <- as.factor("DEAD")
                        CessasionStage <- J
                        CessasionStage <- J
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL                                                                                                                                       
                        
                        if(Fecundity > 0){
                              EPOnset[I] <<- EPStartPing
                              EP[I] <<- Fecundity
                              CEP[I] <<- CBE
                              FirstEprod <- round((((EPStartPing - 1) / (8760 - 1)) - floor(((EPStartPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((EPStartPing - 1) / (8760 - 1)))
                              Termination <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                              FirstEprod_buffer <- FirstEprod - WeightingBuffer
                              Termination_buffer <- Termination + WeightingBuffer
                              if(StartPing >= FirstEprod_buffer & StartPing <= Termination_buffer){ Fitness[I] <<- CurrentFitness }else{ Fitness[I] <<- 0 }
                        }else{
                              Fitness[I] <<- 0
                        }
                        
                        PeakCWeight[I] <<- MaxCWeight
                        EndTime[I] <<- CurrentPing
                        EndStage[I] <<- J
                        break
                  }else if(CPRemap > CessasionHorizon_T){
                        LifeState <- as.factor("KIA")
                        CessasionStage <- J
                        IterationID <- 0                                                                                                                                      
                        DTTrack <- NULL                                                                                                                                       
                        
                        if(Fecundity > 0){
                              EPOnset[I] <<- EPStartPing
                              EP[I] <<- Fecundity
                              CEP[I] <<- CBE
                              FirstEprod <- round((((EPStartPing - 1) / (8760 - 1)) - floor(((EPStartPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((EPStartPing - 1) / (8760 - 1)))
                              Termination <- round((((CurrentPing - 1) / (8760 - 1)) - floor(((CurrentPing - 1) / (8760 - 1)))) * 8760, 0) + 1 - floor(((CurrentPing - 1) / (8760 - 1)))
                              FirstEprod_buffer <- FirstEprod - WeightingBuffer
                              Termination_buffer <- Termination + WeightingBuffer
                              if(StartPing >= FirstEprod_buffer & StartPing <= Termination_buffer){ Fitness[I] <<- CurrentFitness }else{ Fitness[I] <<- 0 }
                        }else{
                              Fitness[I] <<- 0
                        }
                        
                        PeakCWeight[I] <<- MaxCWeight
                        EndTime[I] <<- CurrentPing
                        EndStage[I] <<- J 
                        break
                  }else{
                        #Proceed iteration
                  }
            }
      }
}

