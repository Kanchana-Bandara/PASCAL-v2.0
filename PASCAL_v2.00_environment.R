PhysicalParameters_D <- function(){
      phi <- Latitude
      psi <- WCLAttenCoef
      q <- NVPRFactor
      LIndex <- which(EPSummary[, 1] == phi)
      
      #Irradiance
      jd <- JD(seq(ISOdate(2015, 01, 01, 0), ISOdate(2015, 12, 31, 23), by = "hour"))
      SAngles <- sunpos(sunvector(jd, phi, 15, 1))
      Theta <- as.vector(SAngles[, 2])
      GHI <- NULL
      
      for (i in 1: length(Theta)) {
            Irad <- round((1159.24 * (cos_d(Theta[i]))^1.179) * (exp(-0.0019 * (90 - Theta[i]))), 2)
            GHI <- append(GHI, Irad)
      }
      
      GHI[is.na(GHI)] <- 0
      GHI <- ((GHI - min(GHI)) / (max(GHI) - min(GHI))) * 1000
      
      for(i in 1:8760){
            Irradiance[, i] <<- GHI[i] * exp(psi * MinDepth:MaxDepth)
      }
      
      #Temperature :: SST, SML, TC
      MinTemp <- EPSummary$SSTMIN[LIndex]
      MaxTemp <- EPSummary$SSTMAX[LIndex]
      TSeries <- 1:8760
      mu <- EPSummary$SST_PEAK[LIndex]
      sigma <- EPSummary$SIGMA[LIndex]
      SST <- (1 / (sigma * sqrt(2 * pi))) * (exp(-(((TSeries - mu)^2) / (2 * sigma^2))))
      SST <- (SST - min(SST)) / (max(SST) - min(SST))
      SST <- round(SST * (MaxTemp - MinTemp) + MinTemp, 2)
      
      MLDMin <- 5
      MLDMax <- MaxDepth / 2
      ML_T1 <- EPSummary$SST_S[LIndex] #Time at which ML begins to increase :: sync with SST
      ML_T2 <- EPSummary$PPB_S[LIndex] #Time at which ML becomes shallowest :: sync with PPB onset
      ML_T3 <- EPSummary$PPB_E[LIndex] #Time at which ML begins to get deeper :: sync with PPB termination
      ML_T4 <- EPSummary$SST_E[LIndex] #Time at which ML becomes deepest :: sync with SST
      
      ML_M1 <- rep(MLDMax, ML_T1) #Early deep phase of ML
      MLSRate <- (MLDMin - MLDMax) / (ML_T2 - ML_T1) #Rate of shallowing of ML (m/h)
      ML_M2 <- seq(MLDMax, MLDMin, MLSRate) #Increasing phase of ML
      ML_M3 <- rep(MLDMin, (ML_T3 - ML_T2)) #Stagnenet phase of ML
      MLDRate <- (MLDMax - MLDMin) / (ML_T4 - ML_T3)
      ML_M4 <- seq(MLDMin, MLDMax, MLDRate) #Increasing phase of ML
      ML_M5 <- rep(MLDMax, 8760 - length(c(ML_M1, ML_M2, ML_M3, ML_M4))) #Late deep phase of ML
      MLD <- c(ML_M1, ML_M2, ML_M3, ML_M4, ML_M5) #Complex model
      MLD <- round(MLD, 0)
      MixedLayerDepth <<- MLD
      
      TDRate <- -0.05 #Dissipation rate (°C/m)
      
      for(i in 1:8760){
            SMLTemp <- rep(SST[i], MLD [i])
            TCDTemp <- round(seq(SST[i], MinTemp, TDRate), 2)
            HGZTemp <- rep(MinTemp, (MaxDepth - length(c(SMLTemp, TCDTemp))))
            Temperature[i] <<- c(SMLTemp, TCDTemp, HGZTemp)
      }
      
      #Food availability :: Pelagic and ice-associated production
      MinChla <- 0
      MaxChla <- MaxFA
      PDMax <- 30
      PDepth <- ifelse(Setting == 1, yes = 150, no = 250)
      PDSeries <- seq(0, 1, length.out = PDepth)
      PSStart <- EPSummary$PPB_S[LIndex]
      PSEnd <- EPSummary$PPB_E[LIndex]
      PSPeak <- EPSummary$PPB_PEAK[LIndex]
      PSeries <- seq(0, 1, length.out = (PSEnd - PSStart))
      alpha_PS1 <- 2 
      alpha_PS2 <- 1 / ((PSPeak - PSStart) / (PSEnd - PSStart))
      PS <- (PSeries^(alpha_PS1 - 1) * (1 - PSeries)^(alpha_PS2 - 1) / (beta(alpha_PS1, alpha_PS2)))
      LBank <- rep(0, PSStart)
      RBank <- rep(0, (8760 - PSEnd))
      Chla <- c(LBank, PS, RBank)
      Chla <- (Chla - min(Chla)) / (max(Chla) - min(Chla))
      Chla <- round(Chla * (MaxChla - MinChla) + MinChla, 2)
      
      alpha_PD1 <- 1
      alpha_PD2 <- 3                                                                                  #higher values, more aggrassive depletion
      #alpha_PD2 <- 1 / (PDMax / PDepth) 
      PVD <- (PDSeries^(alpha_PD1 - 1) * (1 - PDSeries)^(alpha_PD2 - 1) / (beta(alpha_PD1, alpha_PD2))) 
      #PVD <- c(PVD[6:110], rep(0, MaxDepth - length(PVD[6:110])))
      PVD <- (PVD - min(PVD)) / (max(PVD) - min(PVD))
      PVD <- c(PVD, rep(0, (MaxDepth - PDepth)))
      
      for(i in 1:8760){
            FA[, i] <<- round(Chla[i] * PVD, 2)
      }
      
      #Predation risk
      IMin <- min(Irradiance)
      IMax <- max(Irradiance)
      
      for(i in 1:8760){
            INorm <- (Irradiance[, i] - IMin) / (IMax - IMin)
            IVPrisk <- INorm * (0.9 - 0.1) + 0.1
            IPrisk <- IVPrisk + q
            PredationRisk[, i] <<- IPrisk
      }
}
