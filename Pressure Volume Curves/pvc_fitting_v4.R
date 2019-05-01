pvc_fitting <- function(total_mass,dry_mass,extra_mass,water_pot,ip,bkp,leafA,LA,ID,data=NULL){
  
  # 1. Assign columns to objects within the function
  
  total_mass <- eval(substitute(total_mass),data, parent.frame())
  dry_mass <- eval(substitute(dry_mass),data, parent.frame())
  extra_mass <- eval(substitute(extra_mass),data,parent.frame())
  water_pot <- eval(substitute(water_pot),data,parent.frame())
  ID <- eval(substitute(ID),data,parent.frame())
  
  # 2. Create matrix to print results
  
  PVC_output <- matrix(nrow = 1,ncol = 13,
                       dimnames = list(c(paste(unique(ID))),
                                       c("IP","MWC.g","R2.MWC","SWC.g","OP.FT.MPa","BP","TLP.MPa","R2.TLP","RWC.TLP.%",
                                         "MEL.MPa","CFT.MPa-1","CTLP.MPa-1","CFT.molm-2MPa-1")))
  
  # 3. Calculate water mass
  
  water_g <- water_mass(total.mass = total_mass,dry.mass = dry_mass,extra.mass = extra_mass)# water mass calculation
  
  # 4. Estimation of 100% water content, use bkp as breakpoint
  
  lm.rwc <- lm(water_g[ip:bkp]~water_pot[ip:bkp])# estimation of 100% Water content
  
  PVC_output[1,1] <- ip
  
  PVC_output[1,2] <- round(lm.rwc$coefficients[1],digits = 3)# 100% Water Content
  
  PVC_output[1,3] <- round(summary(lm.rwc)$r.squared,digits = 2)# 100% water content estimation fit
  
  RWC_p <- ((water_g/lm.rwc$coefficients[1])*100)# Relative water content
  
  RWC <- (water_g/lm.rwc$coefficients[1])
  
  iRWC_p <- 100-RWC_p # Relative water loss
  
  PVC_output[1,4] <- round(PVC_output[1,1]/dry_mass[1],digits = 3)# Saturated water content
  
  # 5. Estimation of Turgor loss point
  
  PVC_output[1,6] <- bkp
  
  lmTLP <- lm((-1/water_pot[bkp:length(iRWC_p)])~iRWC_p[bkp:length(iRWC_p)])# linear fit to estimate TLP
  
  PVC_output[1,7] <- round(-1/(lmTLP$coefficients[1]+lmTLP$coefficients[2]*iRWC_p)[bkp],digits = 3)# Calculate Water potential at TLP
  
  PVC_output[1,8] <- round(summary(lmTLP)$r.squared,digits = 2)# Obtain R2 from TLP estimation
  
  PVC_output[1,9] <- round(RWC_p[bkp],digits = 2)# Relative water content at full turgor
  
  # 6. Estimation of Osmotic potential at full turgor
  
  PVC_output[1,5] <- round(-1/lmTLP$coefficients[1],digits = 3)# Osmotic potential at full turgor
  
  O_WP <- -1/(lmTLP$coefficients[1]+lmTLP$coefficients[2]*iRWC_p)# Osmotic potentials
  
  psiP <- water_pot-O_WP# PSi potential
  
  # 7. Estimation of the modulus of elasticity
  
  lmEl <- lm(psiP[ip:bkp]~(RWC[ip:bkp]))# Estimation of the modulus of elasticity
  
  PVC_output[1,10] <- round(lmEl$coefficients[2],digits = 3)# Modulus of elasticity
  
  # 8. Estimation of the Capacitance at full turgor
  
  lmCFT <- lm(RWC[ip:bkp]~water_pot[ip:bkp])# Linear fit to estimate capacitance at full turgor
  
  PVC_output[1,11] <- round(lmCFT$coefficients[2],digits = 3)# capacitance at full turgor
  
  # 8. Estimation of the Capacitance at turgor loss point
  
  lmCTLP <- lm(RWC[bkp:length(RWC)]~water_pot[bkp:length(RWC)])# Linear fit to estimate capacitance at TLP
  
  PVC_output[1,12] <- round(lmCTLP$coefficients[2],digits = 3)# capacitance at TLP
  
  if(leafA == TRUE){
    PVC_output[1,13] <- round((PVC_output[1,9]*lm.rwc$coefficients[1])/18/(LA/10000),digits = 3)
  }else{
    print("No Leaf Area provided")
  }
  
  par(mfrow=c(1,2))
  
  plot(water_g~water_pot,pch=19,cex=1.2,xlim=c(min(water_pot),0),ylim=c(min(water_g),lm.rwc$coefficients[1]),
       ylab = "Water mass (g)",xlab = "Water Potential (MPa)",
       main=paste(unique(ID),"TLP:",round(-1/(lmTLP$coefficients[1]+lmTLP$coefficients[2]*iRWC_p)[bkp],digits = 2)))
  abline(lm(water_g[ip:bkp]~water_pot[ip:bkp]))
  plot((-1/water_pot)~iRWC_p,pch=19,cex=1.2,
       ylab = "Water Potential (-1/MPa)",xlab = "100-RWC",main=paste("IP:",ip,"(",round(summary(lm.rwc)$r.squared,digits = 2),")",
                                                                     "; BKP:",bkp,"(",round(summary(lmTLP)$r.squared,digits = 2),")"))
  abline(lm((-1/water_pot[bkp:length(iRWC_p)])~iRWC_p[bkp:length(iRWC_p)]))

  
  return(PVC_output)
  
}
