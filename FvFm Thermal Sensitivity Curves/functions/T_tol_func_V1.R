Thermal_tol <- function(fv_fm,temp,species,sp_levels,id,data=NULL,print_graph=TRUE,graph_dir){
  # 1. Assign columns to objects within the function
  
  fv_fm <- eval(substitute(fv_fm),data, parent.frame())
  temp <- eval(substitute(temp),data, parent.frame())
  species <- eval(substitute(species),data,parent.frame())
  id <- eval(substitute(id),data,parent.frame())
  
  data <- data.frame(species,id,fv_fm,temp)# new data frame to perform the analysis
  
  # 2. Create list to print results over all the species and individuals
  
  Thermal_tol_results <- list()
  
  # 3. Loop over species and individuals to run the model
  
  for(i in 1:length(sp_levels)){
    
    datum <- data[data$species == sp_levels[i],]
    
    TT_output <- matrix(nrow = max(datum$id),ncol = 5,
                        dimnames = list(c(paste(unique(datum$id))),
                                        c("ID","Fv_Fm_Max","T50","Tcrit","R2")))
    id_levels <- unique(datum$id)
    
    for(t in 1:max(datum$id)){
      
      datum2 <- datum[datum$id == id_levels[t],]
      
      print(paste(unique(datum2$species),unique(datum2$id),sep="-"))
      
      # 4. Perform linear logit regression to extract initial parameters for theta2 and theta3
      cof <- numeric()
      
      cof[1] <- lm(logit(fv_fm)~temp,data=datum2)$coefficients[1]# theta 2
      
      cof[2] <- lm(logit(fv_fm)~temp,data=datum2)$coefficients[2]# theta 3
      
      # 5. Fit the logistic decay function of FV/Fm as a function of temperature
      
      try(hotmod <- nls(fv_fm ~ theta1/(1 + exp(-(theta2 + theta3*temp))),
                    start= list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
                    data= datum2,
                    trace= F,
                    control = list(maxiter=1000, tol=1e-3)))
      
      # 6. Extract parameters from the non linear model
      
      phia <- coef(hotmod)[1]# fmax, maximum fluorescense
      TT_output[t,2] <- phia
      phib <- coef(hotmod)[2]
      phic <- coef(hotmod)[3]
      
      # 7. Predict Fv/Fm for a range of temperatures based on the model
      x2 <- seq(from= min(datum2$temp),to = max(datum2$temp),by = 0.5)
      
      x3 <- datum2$temp
      
      # 8. Calculate CTmax, or T50
      
      y2<-phia/(1+exp(-(phib+phic*x2))) 
      
      half<- max(y2)*0.50
      
      ctmax.fl1 = (-log((phia/half)-1)-phib)/phic ##***This is CTmax, T50***###
      
      TT_output[t,3] <- ctmax.fl1
      
      # 9. Calculate Tcrit
      
      # 10. Estimate non linear model fit by plotting real Fv/Fm to predicted Fv/Fm
      
      y3 <- phia/(1+exp(-(phib+phic*x3))) 
      
      model.fit <- lm(datum2$fv_fm~y3)
      
      TT_output[t,5] <- summary(model.fit)$r.squared
      
      # 11. Get species and individual name in the output
      
      TT_output[t,1] <- unique(datum2$id)
      
      # 12. Plot the results for a given individual of a given species
      if(!print_graph){
        par(mfrow=c(1,2),oma=c(2,2,2,2), mar=c(2,2,2,2))
        
        plot(datum2$fv_fm~datum2$temp,ylim=c(0,0.88),xlim=c(24,64),type="n")
        text(x = 55,y = 0.88,labels = bquote("F"[v]*"/"*"F"[m]==""))
        text(x = 61,y = 0.88,labels = paste(round(phia,digits = 3)))
        text(x= 55.8,y = 0.82,labels = bquote("T"[50]==""))
        text(x = 61.1,y = 0.82,labels = paste(round(ctmax.fl1,digits = 2)))
        text(x= 55.8,y = 0.76,labels = bquote("T"[crit]==""))
        mtext(text = bquote(Temperature~"("*degree*C*")"),side = 1,line = 2,cex=1)
        mtext(text = bquote("F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
        mtext(text = paste(unique(datum2$species),unique(datum2$id),sep = "-"),side = 3,outer = T)
        lines(x2, y2, col="#B22222",lwd=2.5)
        points(datum2$fv_fm~datum2$temp,pch=16,col="#00000095",cex=1.5)
        
        plot(datum2$fv_fm~y3,ylim=c(0,0.88),xlim=c(0,0.88),type="n")
        text(x = 0.03,y = 0.88,labels = bquote(R^2==""))
        text(x = 0.13,y = 0.872,labels = paste(round(summary(model.fit)$r.squared,digits = 3)))
        abline(a = 0,b = 1,col="#808080",lwd=2.5)
        abline(model.fit,col="#B22222",lwd=2.5)
        mtext(text = bquote(Model~"F"[v]*"/"*"F"[m]),side = 1,line = 2,cex = 1)
        mtext(text = bquote(Data~"F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
        points(datum2$fv_fm~y3,pch=16,col="#00000095",cex=1.5)
      }else{
        jpeg(filename=paste(graph_dir,#directory
                            as.character(unique(datum2$species)),#species_name
                            as.character(unique(datum2$id)),#individual_id
                            ".jpg",sep = ""),width = 5500, height = 2500, units = "px", res = 550)
        par(mfrow=c(1,2),oma=c(2,2,2,2), mar=c(2,2,2,2))
        
        plot(datum2$fv_fm~datum2$temp,ylim=c(0,0.88),xlim=c(24,64),type="n")
        text(x = 55,y = 0.88,labels = bquote("F"[v]*"/"*"F"[m]==""))
        text(x = 61,y = 0.88,labels = paste(round(phia,digits = 3)))
        text(x= 55.8,y = 0.82,labels = bquote("T"[50]==""))
        text(x = 61.1,y = 0.82,labels = paste(round(ctmax.fl1,digits = 2)))
        text(x= 55.8,y = 0.76,labels = bquote("T"[crit]==""))
        mtext(text = bquote(Temperature~"("*degree*C*")"),side = 1,line = 2,cex=1)
        mtext(text = bquote("F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
        mtext(text = paste(unique(datum2$species),unique(datum2$id),sep = "-"),side = 3,outer = T)
        lines(x2, y2, col="#B22222",lwd=2.5)
        points(datum2$fv_fm~datum2$temp,pch=16,col="#00000095",cex=1.5)
        
        plot(datum2$fv_fm~y3,ylim=c(0,0.88),xlim=c(0,0.88),type="n")
        text(x = 0.03,y = 0.88,labels = bquote(R^2==""))
        text(x = 0.13,y = 0.872,labels = paste(round(summary(model.fit)$r.squared,digits = 3)))
        abline(a = 0,b = 1,col="#808080",lwd=2.5)
        abline(model.fit,col="#B22222",lwd=2.5)
        mtext(text = bquote(Model~"F"[v]*"/"*"F"[m]),side = 1,line = 2,cex = 1)
        mtext(text = bquote(Data~"F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
        points(datum2$fv_fm~y3,pch=16,col="#00000095",cex=1.5)
        dev.off()
      }
    }
    
    TT_output <- data.frame(species=unique(datum$species),TT_output)
    
    Thermal_tol_results[[i]] <- TT_output
    
  }
  return(Thermal_tol_results)
}