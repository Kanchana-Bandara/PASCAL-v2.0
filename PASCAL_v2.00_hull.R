#########
#Seeding#
#########

#Dependencies
require(aspace)
require(insol)
require(compiler)

#RNG and paths
SimulationID <- 1                                                                               #Auto update as increments
DirID_1 <- paste0("Simulation_", SimulationID)
DirID_2 <- paste0("D:/Paper_3_inprogress/FinalModel_2/ModelOutputs/", DirID_1)

repeat{
      if(dir.exists(DirID_2)){
            SimulationID <- SimulationID + 1
            DirID_1 <- paste0("Simulation_", SimulationID)
            DirID_2 <- paste0("D:/Paper_3_inprogress/FinalModel_2/ModelOutputs/", DirID_1)
      }else{
            SimulationID <- SimulationID
            DirID_1 <- paste0("Simulation_", SimulationID)
            DirID_2 <- paste0("D:/Paper_3_inprogress/FinalModel_2/ModelOutputs/", DirID_1)
            dir.create(DirID_2)
            set.seed(SimulationID)
            break
      }
}

#Basic input parameters
NSource <- c(2.5e6, 2.5e5)                                                                          #Population sizes of C-0 and C1
DiapauseMachineCode <- as.factor("ON")                                                          #High-accuracy Machine processing of diapause duration dynamics
CZSimulation <- as.factor("ON")                                                                 #Binary C-0 validity
ELimit <- 0.050                                                                                 #How much fraction to be extracted from CZSimulation
N <- ifelse(CZSimulation == "OFF", yes = NSource[2], no = NSource[1])                           #INPUT :: 1-1000000
PBTicks <- 100                                                                                  #No. of ticks in the progressbar
SpeciesID <- 2                                                                                  #INPUT :: 2: CF, 3:CG, 4:CH
CarrierID <- NULL                                                                               #Sequence of individuals as sprawned
Latitude <- 60                                                                                  #INPUT :: 60, 70 or 08
MinDepth <- 1                                                                                   #=1m
Setting <- 2                                                                                    #INPUT :: 1: coastal, 2: open-waters
MaxDepth <- ifelse(Setting == 1, yes = 250, no = 1000)                                          #=250 OR 1000, based on the setting, see above
MaxFA <- 6                                                                                      #INPUT :: > 0
WCLAttenCoef <- -0.06                                                                           #Water column light attenuation coefficient
NVPRFactor <- 0.01                                                                              #Constant non-visual predation risk
AFS <- as.factor("OFF")                                                                         #Binary for alternative food sources
AFSSeasonality <- as.factor("OFF")                                                              #Advanced seasonal distribution of alternative food sources
AFSPeak <- 0.1                                                                                  #Peak alternaticve food source concentration, fraction of MaxFA
AFSBg <- 0.01                                                                                   #Background AFS concentration
AMC_L <- 1                                                                                      #Metabolic cost scalar for movements elsewhere
AMC_H <- 2                                                                                      #Metabolic cost scalar for movements in the surface-mixed layer
DMC <- 0.25                                                                                     #Diapause metabolic cost reduction factor (0.1–0.25)
VPRScalar <- 1e-4                                                                               #Visual predation risk scalar
NVPRScalar <- 0                                                                                 #Non-visual predation risk (already incorporated to the model, hence nullified)
ChlC <- 30                                                                                      #Food quality (Chl-a:C)
StarvationTolerance <- 0.1                                                                      #Scalar for starvation tolerance (fraction of mass lost)
WeightingBuffer <- 24                                                                           #In hours, the time bufferred for fitness weighing

#Physcal environmental sub-model
Temperature <- as.data.frame(matrix(rep(NA, (8760 * MaxDepth)), ncol = 8760))                   #In °C, Kelvin conversion enclosed in Layer-2
Irradiance <- as.data.frame(matrix(rep(NA, (8760 * MaxDepth)), ncol = 8760))                    #umol.m-2.sec-1, standardized to 1000
PredationRisk <- as.data.frame(matrix(rep(NA, (8760 * MaxDepth)), ncol = 8760))                 #Probabilities in the range of 0.91-0.11
FA <- as.data.frame(matrix(rep(NA, (8760 * MaxDepth)), ncol = 8760))                            #mgChl-a.m-3, ugC.l-1 conversion enclosed
MixedLayerDepth <- NULL                                                                         #For swimming cost estimations, see above
EPSummary <- read.table("EPSummary.txt", sep = " ")

PhysicalParameters_D()

#Limits and constraints
NDim <- 6                                                                                       #No. of dimensions to optimize
STL <- EPSummary$PPB_PEAK[EPSummary$LAT == Latitude] - (90 * 24)                                #Lower-bound for STime :: 90-d prior to the pps onset
STL <- ifelse(STL <= 1, yes = 1, no = STL)                                                      #Lower bound correction
STU <- EPSummary$PPB_PEAK[EPSummary$LAT == Latitude] + (120 * 24)                               #Upper-bound for STime :: 120-d after the pps onset
STU <- ifelse(STU >= 8760, yes = 8760, no = STU)                                                #Upper bound correction
ISPL <- 1e-4                                                                                    #Lower-bound for ISP :: 1e-4 umol.m-2.s1
ISPU <- 1000                                                                                    #Lower-bound for ISP :: 1000 umol.m-2.s1
GAPL <- 0.1                                                                                     #Lower-bound for GAP
GAPU <- 1.0                                                                                     #Upper-bound for GAP
SDPL <- 0.1                                                                                     #Lower-bound for SDP
SDPU <- 0.7                                                                                     #Lower-bound for SDP
SAPL <- 0.1                                                                                     #Lower-bound for SAP
SAPU <- 1.0                                                                                     #Lower-bound for SAP
LCLPL <- 0                                                                                      #Lower-bound for open parameter :: free for prospective use
LCLPU <- 1                                                                                      #Lower-bound for open parameter :: free for prospective use

#Body mass limit implementation
SSGDParameters <- read.table("SSGDParameters.txt", header = TRUE)
CM_Upper <- NULL                                                                                #Vector of maximum critical molting masses, J = 0 to 12
CM_Lower <- NULL                                                                                #Vector of minimum critical molting masses, J = 0 to 12

CritMoltMass_UPR()                                                                              #Function of robust UPPER molting mass formulations :: min temperature and max.food concentration
CritMoltMass_LWR()                                                                              #Function of robust LOWER molting mass formulations :: max temperature and max.food concentration

#Breaking condition
LCBreak <- EPSummary$PPB_E[which(EPSummary$LAT == Latitude)]                                    #Life cycle breaking threshold, in terms of time units :: at the end of PPS

#Baseline vector for allometric scaling of ISP
AS_x <- seq(0, 2800, 1)
AS_a <- 2
AS_b <- 600
AS_c <- 400
AS_y <- AS_a / (1 + exp((AS_b - AS_x)/AS_c))
AS_ynorm <- (AS_y - min(AS_y)) / (max(AS_y) - min(AS_y))
AS_yext <- AS_ynorm * (2500 - 1) + 1
ISP_base <- AS_yext
ISP_base_rc <- 1 / ISP_base

#Baseline vector for allometric scaling of K
AS_x <- seq(0, 2800, 1)
AS_a <- 2
AS_b <- 600
AS_c <- 400
AS_y <- AS_a / (1 + exp((AS_b - AS_x)/AS_c))
AS_ynorm <- (AS_y - min(AS_y)) / (max(AS_y) - min(AS_y))
AS_yext <- AS_ynorm * (2.5 - 0) + 0
K_base <- AS_yext

rm(AS_x, AS_y, AS_a, AS_b, AS_c, AS_ynorm, AS_yext)

#Complex function for storage volume scaling with body mass :: only for dynamic body masses
SV_x <- seq(38, 159, 1)       
SV_a <- 1
SV_b <- 60
SV_c <- 20
SV_2 <- SV_a / (1 + exp((SV_b - SV_x)/SV_c))
SV_2_norm <- (SV_2 - min(SV_2)) / (max(SV_2) - min(SV_2))
SV_1 <- rep(0, 38)
SV_3 <- rep(1, (2800 - 160))
SVol <- round(c(SV_1, SV_2_norm, SV_3), 2)
rm(SV_x, SV_b, SV_a, SV_c, SV_2, SV_2_norm, SV_1, SV_3)

#Condition-zero simulation
if(CZSimulation == "OFF"){
      #Parameter pools
      STimePool <- seq(STL, STU, 24)
      ISPPool <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
      GAPPool <- seq(GAPL, GAPU, 0.1)                                                                 #Modified for LX-MPTM, precision is built internally
      SDPPool <- seq(SDPL, SDPU, 0.1)                                                                 #Modified for LX-MPTM, precision is built internally
      SAPPool <- seq(SAPL, SAPU, 0.1)                                                                 #Modified for LX-MPTM, precision is built internally
      LCLPPool <- seq(LCLPL, LCLPU, 0.1)
      
      #Evolvable (soft) parameters (6D)
      STime <- sample(STimePool, N, replace = TRUE)                                                   #Time of birth, 24H intervals
      isp <- sample(ISPPool, N, replace = TRUE)                                                       #Irradiance Sensitivity Parameter, x10 intervals
      gap <- sample(GAPPool, N, replace = TRUE)                                                       #Growth Allocation Parameter, 0.01 fine intervals
      sdp <- sample(SDPPool, N, replace = TRUE)                                                       #Seasonal Descent Parameter, 0.01 fine intervals
      sap <- sample(SAPPool, N, replace = TRUE)                                                       #Seasonal Ascent Parameter, 0.01 fine intervals
      lclp <- sample(LCLPPool, N, replace = TRUE)                                                     #Life Cycle Length Parameter, binary (for now)
      
      #Proceed to C-1 directly
      
      #Output loggers
      DEntry <- rep(NA, N)
      DExit <- rep(NA, N)
      DStage <- rep(NA, N)
      DDepth <- rep(NA, N)
      DCWeight <- rep(NA, N)
      DEntrySWeight <- rep(NA, N)
      DExitSWeight <- rep(NA, N)
      EPOnset <- rep(NA, N)
      EP <- rep(NA, N)
      CEP <- rep(NA, N)
      Fitness <- rep(NA, N)
      PeakCWeight <- rep(NA, N)
      EndTime <- rep(NA, N)
      EndStage <- rep(NA, N)
      
      #Bytecode compilation
      Simulator_1YGT_c <- cmpfun(Simulator_1YGT_3GV1)
      
      ####################
      #Initial simulation#
      ####################
      ST <- Sys.time()
      PB <- txtProgressBar(min = 0, max = N, style = 3, width = PBTicks)
      BreakPoints <- c(1, seq(100, N, by = 100))
      
      for(CarrierID in 1:N){
            Simulator_1YGT_c()
            
            if(any(CarrierID == BreakPoints)){
                  Sys.sleep(0.1)
                  setTxtProgressBar(PB, CarrierID)
            }else{
                  #proceed
            }
      }
      
      ET <- Sys.time()
      print(ET - ST)
      Fitness[is.na(Fitness)] <- 0 
      Output <- rbind(STime, isp, gap, sdp, sap, lclp, DEntry, DExit, DStage, DDepth, DCWeight, DEntrySWeight, DExitSWeight, EPOnset, EP, CEP, PeakCWeight, EndStage, EndTime, Fitness)
      Output <- t(Output)
      write.table(Output, file = paste0(DirID_2, "/Generation_0.txt"), sep = " ", row.names = FALSE, col.names = TRUE)
}else{
      #Parameter pools
      STimePool <- seq(STL, STU, 24)
      ISPPool <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
      GAPPool <- seq(GAPL, GAPU, 0.01)                                                                 #Modified for LX-MPTM, precision is built internally
      SDPPool <- seq(SDPL, SDPU, 0.01)                                                                 #Modified for LX-MPTM, precision is built internally
      SAPPool <- seq(SAPL, SAPU, 0.01)                                                                 #Modified for LX-MPTM, precision is built internally
      LCLPPool <- seq(LCLPL, LCLPU, 0.1)
      
      #Evolvable (soft) parameters (6D)
      STime <- sample(STimePool, N, replace = TRUE)                                                   #Time of birth, 24H intervals
      isp <- sample(ISPPool, N, replace = TRUE)                                                       #Irradiance Sensitivity Parameter, x10 intervals
      gap <- sample(GAPPool, N, replace = TRUE)                                                       #Growth Allocation Parameter, 0.01 fine intervals
      sdp <- sample(SDPPool, N, replace = TRUE)                                                       #Seasonal Descent Parameter, 0.01 fine intervals
      sap <- sample(SAPPool, N, replace = TRUE)                                                       #Seasonal Ascent Parameter, 0.01 fine intervals
      lclp <- sample(LCLPPool, N, replace = TRUE)                                                     #Life Cycle Length Parameter, binary (for now)
      
      #Proceed to C-0
      
      #Output loggers
      DEntry <- rep(NA, N)
      DExit <- rep(NA, N)
      DStage <- rep(NA, N)
      DDepth <- rep(NA, N)
      DCWeight <- rep(NA, N)
      DEntrySWeight <- rep(NA, N)
      DExitSWeight <- rep(NA, N)
      EPOnset <- rep(NA, N)
      EP <- rep(NA, N)
      CEP <- rep(NA, N)
      Fitness <- rep(NA, N)
      PeakCWeight <- rep(NA, N)
      EndTime <- rep(NA, N)
      EndStage <- rep(NA, N)
      
      #Bytecode compilation
      Simulator_1YGT_c <- cmpfun(Simulator_1YGT_3GV1)
      
      ####################
      #Initial simulation#
      ####################
      ST <- Sys.time()
      PB <- txtProgressBar(min = 0, max = N, style = 3, width = PBTicks)
      BreakPoints <- c(1, seq(100, N, by = 100))
      
      for(CarrierID in 1:N){
            Simulator_1YGT_c()
            
            if(any(CarrierID == BreakPoints)){
                  Sys.sleep(0.1)
                  setTxtProgressBar(PB, CarrierID)
            }else{
                  #proceed
            }
      }
      
      ET <- Sys.time()
      print(ET - ST)
      Fitness[is.na(Fitness)] <- 0 
      PreOutput <- rbind(STime, isp, gap, sdp, sap, lclp, DEntry, DExit, DStage, DDepth, DCWeight, DEntrySWeight, DExitSWeight, EPOnset, EP, CEP, PeakCWeight, EndStage, EndTime, Fitness)
      PreOutput <- as.data.frame(t(PreOutput))
      write.table(PreOutput, file = paste0(DirID_2, "/Generation_0.txt"), sep = " ", row.names = FALSE, col.names = TRUE)
      
      N <- NSource[2]
      PreOutput <- PreOutput[order(PreOutput$Fitness, decreasing = TRUE), ]
      
      if(all(PreOutput$Fitness == 0)){
            Sys.sleep(0.1)
            writeLines("\n***EXECUTION HALTED :: All strategies posess zero fitness***\n")
            Sys.sleep(10^6)
      }else{
            #Basic selection (high-pressure selection)
            ExtractionLimit <- N * ELimit
            MaxPosFitness <- length(which(PreOutput$Fitness > 0))
            
            if(MaxPosFitness >= ExtractionLimit){
                  Extraction_1_Set <- sample(which(PreOutput$Fitness > 0), ExtractionLimit, replace = FALSE)
                  Extraction_1 <- as.data.frame(PreOutput[Extraction_1_Set, ])
            }else{
                  Extraction_1 <- as.data.frame(PreOutput[1:ExtractionLimit, ])
            }
            
            Extraction_2_Set <- sample(ExtractionLimit:nrow(PreOutput), (N - ExtractionLimit), replace = FALSE)
            Extraction_2 <- as.data.frame(PreOutput[Extraction_2_Set, ])
            Output <- rbind(Extraction_1, Extraction_2)
            RandomSortOrder <- sample(1:nrow(Output), nrow(Output), replace = FALSE)
            Output <- Output[RandomSortOrder, ]
      }
}

########################
#Optimization Algorithm#
########################

#Summarizing
NGen <- 0
TGen <- ifelse(CZSimulation == "ON", yes = 150, no = 100)                                                                                 #50% space allocation for buffering
T1Size <- 2
T2Size <- 10 
CX <- 0.3                                                                                                                                 #Crossover prob. = 1 - CX
LP_a <- 0
LP_b <- 0.5
Pool <- c("STimePool", "ISPPool", "GAPPool", "SDPPool", "SAPPool", "LCLPPool")
MRate <- 0.20
mp <- 2
MDPres <- 1                                                                                                                               #Decimal place precision for MRates < 0.1
SLog <- matrix(rep(NA, (TGen * ncol(Output))), ncol = ncol(Output))
rownames(SLog) <- 1:TGen
colnames(SLog) <- colnames(Output)

repeat{
      NGen <- NGen + 1
      FName <- paste0("/Generation_", NGen, ".txt")
      write.table(Output, file = paste0(DirID_2, FName), sep = " ", row.names = FALSE, col.names = TRUE)
      SLogOut <- apply(Output, 2, mean, na.rm = TRUE)
      SLogOut[c(1, 2, 6, 7, 8, 9, 10, 14, 15, 16, 18, 19)] <- round(SLogOut[c(1, 2, 6, 7, 8, 9, 10, 14, 15, 16, 18, 19)], 0)
      SLogOut[c(3, 4, 5, 11, 12, 13, 17, 20)] <- round(SLogOut[c(3, 4, 5, 11, 12, 13, 17, 20)], 2)
      SLog[NGen, ] <- as.vector(SLogOut)
      
      Sys.sleep(0.1)
      writeLines(sprintf("\nGeneration no: %d\n######################", NGen))
      writeLines(sprintf("Mean Fitness: %f", mean(Output[, ncol(Output)], na.rm = TRUE)))
      writeLines(sprintf("Max. Fitness: %f", max(Output[, ncol(Output)], na.rm = TRUE)))
      writeLines(sprintf("Min. Fitness: %f", min(Output[, ncol(Output)], na.rm = TRUE)))
      writeLines("\nOutput Summary (Most recent 5 generations) : \n")
      
      if(NGen <= 5){
            print(SLog[1:NGen, ])
      }else{
            print(SLog[(NGen - 4):NGen, ])
      }
      
      writeLines("\n\n")
      
      #Parent selection :: Deterministic tournament
      PreGenOtput <- Output
      Population_S0 <- as.matrix(Output[, c(1:NDim, ncol(Output))])                                                                       #Selects the evolvable set of parameters and fitness from the rest
      Population_S1 <- matrix(rep(NA, (ncol(Population_S0) * nrow(Population_S0))), nrow = nrow(Population_S0))                           #Output matrix for the selected parents
      rm(Output)
      gc()
      
      PSCounter <- 0
      
      repeat{
            PSCounter <- PSCounter + 1
            ISelect <- sample(1:nrow(Population_S0), T1Size, replace = TRUE)
            MaxFitRange <- which(Population_S0[ISelect, ncol(Population_S0)] == max(Population_S0[ISelect, ncol(Population_S0)]))
            
            if(length(MaxFitRange) > 1){
                  MaxFitInd <- ISelect[MaxFitRange[round(runif(1, min = 1, max = length(MaxFitRange)), 0)]]
                  Population_S1[PSCounter, ] <- Population_S0[MaxFitInd, ]
            }else{
                  MaxFitInd <- ISelect[MaxFitRange]
                  Population_S1[PSCounter, ] <- Population_S0[MaxFitInd, ]
            }
            
            if(PSCounter == (nrow(Population_S0))){
                  break
            }else{
                  #Proceed selection
            }
      }
      
      #Random mating
      Population_S1 <- Population_S1[, -ncol(Population_S1)] #Fitness column is dropped
      PaRange <- 1:nrow(Population_S1)
      FaRange <- sample(1:length(PaRange), length(PaRange)/ 2, replace = FALSE)
      MaRange <- PaRange[-c(FaRange)]
      
      #La-place crossover (LX)
      Population_S2_Fa <- matrix(rep(NA, (length(FaRange) * ncol(Population_S1))), ncol = ncol(Population_S1))
      Population_S2_Ma <- matrix(rep(NA, (length(MaRange) * ncol(Population_S1))), ncol = ncol(Population_S1))
      
      for(i in 1:length(FaRange)){
            Fa <- Population_S1[FaRange[i], ]
            Ma <- Population_S1[MaRange[i], ]
            
            for(j in 1:ncol(Population_S1)){
                  CXProb <- round(runif(1), 1)
                  
                  if(CXProb >= CX){
                        #Do crossover following the LX algorithm
                        LPRand <- round(runif(1, min = 0.1, max = 0.9), 1)
                        if(LPRand <= 0.5){ LPPar <- LP_a - (LP_b * log(LPRand)) }else{ LPPar <- LP_a  + (LP_b * log(LPRand)) }
                        Da1 <- Fa[j] + (LPPar * abs(Fa[j] - Ma[j]))
                        CompPool <- eval(parse(text = Pool[j]))
                        
                        if(Da1 < min(CompPool)){
                              Da1 <- min(CompPool)
                        }else if(Da1 > max(CompPool)){
                              Da1 <- max(CompPool)
                        }else{
                              Da1 <- Da1
                        }
                        Population_S2_Fa[i, j] <- Da1
                        
                        Da2 <- Ma[j] + (LPPar * abs(Fa[j] - Ma[j]))
                        CompPool <- eval(parse(text = Pool[j]))
                        
                        if(Da2 < min(CompPool)){
                              Da2 <- min(CompPool)
                        }else if(Da2 > max(CompPool)){
                              Da2 <- max(CompPool)
                        }else{
                              Da2 <- Da2
                        }
                        Population_S2_Ma[i, j] <- Da2
                  }else{
                        #Proceed with no crossover
                        Da1 <- Fa[j]
                        Population_S2_Fa[i, j] <- Da1
                        Da2 <- Ma[j]
                        Population_S2_Ma[i, j] <- Da2
                  }
            }
      }
      
      Population_S2 <- rbind(Population_S2_Fa, Population_S2_Ma)
      Population_S2[, c(1, 2)] <- round(Population_S2[, c(1, 2)], 0)                                                                            #Zero precision for STime and isp
      Population_S2[, 6] <- round(Population_S2[, 6], 1)                                                                                        #Single precision for STime and isp
      Population_S2[, c(3, 4, 5)] <- round(Population_S2[, c(3, 4, 5)], 2)                                                                      #Dual-decimal precision for gap, sdp and sap
      
      rm(Population_S2_Ma, Population_S2_Fa, Population_S1)
      gc()
      
      #Mutation
      for(i in 1:nrow(Population_S2)){
            MProb <- round(runif(ncol(Population_S2)), MDPres)
            
            for(j in 1:ncol(Population_S2)){
                  if(MProb[j] <= MRate){
                        #Uniform mutation (random replacement)
                        CompPool <- eval(parse(text = Pool[j]))
                        Mutant <- sample(CompPool, 1, replace = TRUE)
                        Population_S2[i, j] <- Mutant
                  }else{
                        #No mutation
                  }
            }
      }
      
      Population_S2[, c(1, 2)] <- round(Population_S2[, c(1, 2)], 0)                                                                            #Zero precision for STime and isp
      Population_S2[, 6] <- round(Population_S2[, 6], 1)                                                                                        #Single precision for STime and isp
      Population_S2[, c(3, 4, 5)] <- round(Population_S2[, c(3, 4, 5)], 2)                                                                      #Dual-decimal precision for gap, sdp and sap
      
      #Secondary life cycle simulation
      N <- nrow(Population_S2)
      STime <- Population_S2[, 1]
      isp <- Population_S2[, 2]
      gap <- Population_S2[, 3]
      sdp <- Population_S2[, 4]
      sap <- Population_S2[, 5]
      lclp <- Population_S2[, 6]
      
      DEntry <- rep(NA, N)
      DExit <- rep(NA, N)
      DStage <- rep(NA, N)
      DDepth <- rep(NA, N)
      DCWeight <- rep(NA, N)
      DEntrySWeight <- rep(NA, N)
      DExitSWeight <- rep(NA, N)
      EPOnset <- rep(NA, N)
      EP <- rep(NA, N)
      CEP <- rep(NA, N)
      Fitness <- rep(NA, N)
      PeakCWeight <- rep(NA, N)
      EndTime <- rep(NA, N)
      EndStage <- rep(NA, N)
      
      ST <- Sys.time()
      PB <- txtProgressBar(min = 0, max = N, style = 3, width = PBTicks)
      BreakPoints <- c(1, seq(100, N, by = 100))
      
      for(CarrierID in 1:N){
            Simulator_1YGT_c()
            
            if(any(CarrierID == BreakPoints)){
                  Sys.sleep(0.1)
                  setTxtProgressBar(PB, CarrierID)
            }else{
                  #proceed
            }
      }
      
      ET <- Sys.time()
      print(ET - ST)
      Fitness[is.na(Fitness)] <- 0
      Output <- rbind(STime, isp, gap, sdp, sap, lclp, DEntry, DExit, DStage, DDepth, DCWeight, DEntrySWeight, DExitSWeight, EPOnset, EP, CEP, PeakCWeight, EndStage, EndTime, Fitness)
      Output <- t(Output)
      
      #Suvivor selection :: Round-Robin Tournament
      Population_S3 <- rbind(PreGenOtput, Output)
      rm(Output)
      gc()
      
      RRTScore <- rep(NA, nrow(Population_S3))
      
      for(i in 1:nrow(Population_S3)){
            LocalScore <- rep(NA, T2Size)
            RRTOpponents <- sample(1:nrow(Population_S3), T2Size, replace = FALSE)
            
            for(j in 1:length(RRTOpponents)){
                  if(Population_S3[i, ncol(Population_S3)] > Population_S3[RRTOpponents[j], ncol(Population_S3)]){
                        LocalScore[j] <- 1
                  }else if(Population_S3[i, ncol(Population_S3)] == Population_S3[RRTOpponents[j], ncol(Population_S3)]){
                        LocalScore[j] <- 0
                  }else{
                        LocalScore[j] <- -1
                  }
            }
            
            RRTScore[i] <- sum(LocalScore)
      }
      
      Population_S3 <- cbind(Population_S3, RRTScore)
      Population_S3 <- Population_S3[order(Population_S3[, ncol(Population_S3)], decreasing = TRUE), ]
      Output <- Population_S3[1:(nrow(Population_S3) / 2), -ncol(Population_S3)]
      rm(Population_S2, Population_S3, RRTOpponents, LocalScore, RRTScore)
      gc()
      
      if(NGen >= TGen){
            FName <- paste0("/OptimizedOutput_SIM", SimulationID, ".txt")
            write.table(Output, file = paste0(DirID_2, FName), sep = " ", row.names = FALSE, col.names = TRUE)
            FName <- paste0("/SummaryLog_SIM", SimulationID, ".txt")
            write.table(SLog, file = paste0(DirID_2, FName), sep = " ", row.names = FALSE, col.names = TRUE)
            break
      }
}