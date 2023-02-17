CritMoltMass_LWR <- function(){
      SSGDParameters_Mapsetal_2011 <- read.table("SSGDParameters_mapsetal_2011.txt", header = TRUE)
      CurrentTemp <- max(as.matrix(Temperature)) + 273.15
      T0 <- 273
      k <- 8.62e-5
      a <- 0.6
      S <- SpeciesID
      J <- 3
      
      E_M <- SSGDParameters_Mapsetal_2011[2, S]
      E_D <- SSGDParameters_Mapsetal_2011[4, S]
      E_I <- SSGDParameters_Mapsetal_2011[5, S]
      M0 <- SSGDParameters_Mapsetal_2011[3, S]
      CurrentCWeight <- SSGDParameters_Mapsetal_2011[1, S]
      CurrentFA <- max(FA) * 30
      Counter <- 0
      DTCounter <- 0
      DTTrack <- NULL
      BMass <- NULL
      FLim <- NULL
      JDT <- NULL
      JT <- NULL
      BMassMax <- rep(NA, 13)
      
      repeat{
            Counter <- Counter + 1
            DTCounter <- DTCounter + 1
            JS <- ifelse(J <= 6, yes = 1, no = 2)
            H0S <- SSGDParameters_Mapsetal_2011[6:7, S]
            V0S <- SSGDParameters_Mapsetal_2011[8:9, S]
            F0S <- SSGDParameters_Mapsetal_2011[10:11, S]
            H0 <- H0S[JS]
            V0 <- V0S[JS]
            F0 <- F0S[JS]
            CurrentHandlingTime <- H0 * CurrentCWeight^-0.75 * exp((-E_I * (CurrentTemp - T0)) / (k * CurrentTemp * T0))
            CurrentClearenceRate <- V0 * CurrentCWeight^0.75 * exp((E_I * (CurrentTemp - T0)) / (k * CurrentTemp * T0))
            CurrentIngestionRate <- (CurrentFA * CurrentClearenceRate) / (1 + (CurrentHandlingTime * CurrentFA * CurrentClearenceRate))
            CurrentFL <- CurrentIngestionRate * CurrentHandlingTime
            CurrentIngestionRate <- CurrentIngestionRate * 3.6e-3
            CurrentBMRate <- ((M0 * CurrentCWeight^0.75 * exp((E_M * (CurrentTemp - T0)) / (k * CurrentTemp * T0)))) * 60 * 60
            CurrentAMRate <- 0
            CurrentGPotential <- ((a * CurrentIngestionRate) - (CurrentBMRate + CurrentAMRate))                                                                                     
            CurrentCWeight <- CurrentCWeight + CurrentGPotential
            D0 <- SSGDParameters_Mapsetal_2011[(J + 12), S]
            CurrentDT <- 1 / ((1 / (round(D0 * exp((-E_D * (CurrentTemp - T0)) / (k * CurrentTemp * T0)) * 24, 2))) * CurrentFL^F0)
            
            BMass <- append(BMass, CurrentCWeight)
            FLim <- append(FLim, CurrentFL)
            DTTrack <- append(DTTrack, CurrentDT)
            
            if(DTCounter >= gm_mean(DTTrack)){
                  J <- J + 1
                  JDT <- append(JDT, DTCounter)
                  JT <- append(JT, Counter)
                  DTCounter <- 0
                  DTTrack <- NULL
                  if(J == 12){
                        break
                  }
            }else{
                  #Re-iterate
            }
      }
      
      #Exporting to .GlobalEnv
      BMassLS <- round(BMass[JT], 2)
      BMassMax <- c(rep(SSGDParameters_Mapsetal_2011[1, S], 3), BMassLS, BMassLS[length(BMassLS)])
      CM_Lower <<- BMassMax 
}