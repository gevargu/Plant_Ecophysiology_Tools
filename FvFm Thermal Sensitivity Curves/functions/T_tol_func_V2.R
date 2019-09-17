Thermal_tol <- function(fv_fm,temperature,species,sp_levels,individual,data=NULL,print_graph=TRUE,graph_dir){
  
  require(car)
  
  # 1. Assign columns to objects within the function
  
  fv_fm <- eval(substitute(fv_fm),data, parent.frame())
  temperature <- eval(substitute(temperature),data, parent.frame())
  species <- eval(substitute(species),data,parent.frame())
  individual <- eval(substitute(individual),data,parent.frame())
  
  data <- data.frame(species,individual,fv_fm,temperature)# new data frame to perform the analysis
  
  # 2. Create list to print results over all the species and individuals
  
  Thermal_tol_results <- list()
  
  # 3. Loop over species and individuals to run the model
  
  for(i in 1:length(sp_levels)){
    
    datum <- data[data$species == sp_levels[i],]
    
    ID=species=Fv_Fm_Max=ctmax=ctmax.lci=ctmax.uci=tcrit=tcrit.lci=tcrit.uci=R2=c()
    
    id_levels <- unique(datum$individual)
    
    for(t in 1:max(datum$individual)){
      
      datum2 <- datum[datum$individual == id_levels[t],]
      
      IT <- min(datum2$temperature)
      
      if(datum2[which(datum2$temperature==IT),3] < 0.7){
        
        # 4. Skipping poor quality data, plants with fv/fm < 0.7
        
        print(paste("Skipping ",unique(datum2$species),"-",unique(datum2$individual),": initial fv/fm < 0.7",sep=""))
        
        }else{
          
          # 5. Perform linear logit regression to extract initial parameters for theta2 and theta3
          
          cof <- numeric()
          
          cof[1] <- lm(logit(fv_fm)~temperature,data=datum2)$coefficients[1]# theta 2
          cof[2] <- lm(logit(fv_fm)~temperature,data=datum2)$coefficients[2]# theta 3
          
          # 6. Fit the logistic decay function of FV/Fm as a function of temperature
          
          try(hotmod <- nls(fv_fm ~ theta1/(1 + exp(-(theta2 + theta3*temperature))),
                          start= list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
                          data= datum2,
                          trace= F,
                          control = list(maxiter=1000, tol=1e-3)))
          
          # 7. Extract parameters from the non linear model
          
          phia <- coef(hotmod)[1]# fmax, maximum fluorescense
          Fv_Fm_Max[t] <- phia
          phib <- coef(hotmod)[2]
          phic <- coef(hotmod)[3]
          
          x <- seq(from= min(datum2$temperature),to = max(datum2$temperature),by = 0.5)
          lx <- length(x)
          y <- phia/(1+exp(-(phib+phic*x))) 
          predict <- data.frame(x,y)
          
          #calculates half of the initial Fv/Fm
          initial <- subset(datum2, datum2$temperature==IT, select=fv_fm)
          half <- mean(initial[,1])/2  
          
          # 8. Estimate parameters with 95 Confidence Interval
          predict.boot <- matrix(NA,41,500)
          ctmax.half.bs=Tcrit.bs=c()
          ptm=proc.time()
          
          for(k in 1:100){#start loop for bootstrap estimation
            srows <- sample(1:nrow(datum2),nrow(datum2),TRUE)
            
            cof2 <- coef(lm(logit(fv_fm)~temperature,data=datum2[srows,])) 
            
            if(class(try(nls(fv_fm ~ theta1/(1 + exp(-(theta2 + theta3*temperature))),
                             start=list(theta1 = .8, theta2 = cof2[1], theta3 = cof2[2]),
                             data=datum2[srows,],
                             trace=F,
                             control=list(maxiter=1000, tol=1e-3)), silent=T)[[1]])=="nlsModel"){
              hotmod2 <- nls(fv_fm ~ theta1/(1 + exp(-(theta2 + theta3*temperature))),
                            start=list(theta1 = .8, theta2 = cof2[1], theta3 = cof2[2]),
                            data=datum2[srows,],
                            trace=F,
                            control=list(maxiter=1000, tol=1e-3))
              phia2 <- coef(hotmod2)[1]
              phia2 <- coef(hotmod2)[1]
              phib2 <- coef(hotmod2)[2]
              phic2 <- coef(hotmod2)[3]
              
              x2 <- seq(from= min(datum2$temperature),to = max(datum2$temperature))
              y2 <- phia2/(1+exp(-(phib2+phic2*x2)))
              
              predict.boot[,k] <- y2
              ctmax.fl1 <- (-log((phia2/half)-1)-phib2)/phic2 #estimate ctmax
              ctmax.half.bs[k] <- ctmax.fl1 #estimate ctmax
              predict <- data.frame(x2,y2) #create the prediction data frame
              df1 <- cbind(predict[-1,], predict[-nrow(predict),])
              df1 <- df1[,c(3,1,4,2)]
              df1$slp <- as.vector(apply(df1, 1, function(x) summary(lm((x[3:4]) ~ x[1:2])) [[4]][[2]] ))
              max.slp <- round(min(df1$slp), 3)
              slp.at.tcrit <- max.slp*.15
              b4ct <- which(df1[,1] < ctmax.fl1)
              fvfv.at.tcrit <- df1[which(abs(df1[b4ct,]$slp-slp.at.tcrit)==min(abs(df1[b4ct,]$slp-slp.at.tcrit))),][1,3]
              Tcrit.bs[k] <- round((-log((phia/fvfv.at.tcrit)-1)-phib)/phic, 2)
              }else{
                (class(try(nls(fv_fm ~ theta1/(1 + exp(-(theta2 + theta3*temperature))),
                               start=list(theta1 = .8, theta2 = cof2[1], theta3 = cof2[2]),
                               data=datum2[srows,],
                               trace=F,
                               control=list(maxiter=1000, tol=1e-3)), silent=T)[[1]])=="list")
                predict.boot[,k]="NA"
                ctmax.half.bs[k]="NA"
                Tcrit.bs[k]="NA"
              }
            }# finish loop for bootstrap estimation
          
          # print time ellapsed in bootstrapping
          print(proc.time() - ptm); print(paste(unique(datum2$species),"-",unique(datum2$individual),
                                                ": initial fv/fm > 0.7",sep=""))
          
          predict.boot <- apply(predict.boot,2, as.numeric)
          ctmax.half.bs <- as.numeric(ctmax.half.bs)
          Tcrit.bs <- as.numeric(Tcrit.bs)
          
          #warnings()
          
          #predict.CI=matrix(NA,lx,2)
          
          #for(k in 1:lx){
          #  predict.CI[k,]=quantile(predict.boot[k,],c(0.025,0.975),na.rm=T)
           # }
          
          #ctmax.CI=quantile(ctmax.bs,c(0.025,0.975),na.rm=T)
          ctmax.half.CI <- quantile(ctmax.half.bs,c(0.025,0.975),na.rm=T)
          Tcrit.CI <- quantile(Tcrit.bs,c(0.025,0.975),na.rm=T)
          #Tcrit.CI.slp=quantile(Tcrit.bs,c(0.025,0.975),na.rm=T)
          
          species[t] <- as.character(unique(datum2$species))
          ID[t] <- as.character(unique(datum2$individual))# get individual name
          ctmax[t] <- round(mean(na.omit(ctmax.half.bs)),2)  #predicted cTmax at half of initial Fv/Fm
          ctmax.lci[t] <- round(ctmax.half.CI[1],2)#Ctmax lower confidence interval
          ctmax.uci[t] <- round(ctmax.half.CI[2],2)#Ctmax upper confidence interval
          
          tcrit[t] <- round(mean(na.omit(Tcrit.bs)),2)# predicted tcrit
          tcrit.lci[t]=round(na.omit(Tcrit.CI[1]),2)#tcrit lower confidence interval
          tcrit.uci[t]=round(na.omit(Tcrit.CI[2]),2)#tcrit upper confidence interval

          # 10. Estimate non linear model fit by plotting real Fv/Fm to predicted Fv/Fm
          x3 <- datum2$temperature
          
          y3 <- phia/(1+exp(-(phib+phic*x3))) 
          
          model.fit <- lm(datum2$fv_fm~y3)
          
          R2[t] <- summary(model.fit)$r.squared
          
          # 12. Plot the results for a given individual of a given species
          if(!print_graph){# start plotting action 1
            par(mfrow=c(1,2),oma=c(2,2,2,2), mar=c(2,2,2,2))
            
            plot(datum2$fv_fm~datum2$temperature,ylim=c(0,0.88),xlim=c(24,64),type="n")
            text(x = 55,y = 0.88,labels = bquote("F"[v]*"/"*"F"[m]==""))
            text(x = 61,y = 0.88,labels = paste(round(phia,digits = 3)))
            text(x= 55.8,y = 0.82,labels = bquote("T"[50]==""))
            text(x = 61.1,y = 0.82,labels = paste(round(ctmax[t],digits = 2)))
            text(x= 55.8,y = 0.76,labels = bquote("T"[crit]==""))
            text(x = 61.1,y = 0.76,labels = paste(round(tcrit[t],digits = 2)))
            mtext(text = bquote(Temperature~"("*degree*C*")"),side = 1,line = 2,cex=1)
            mtext(text = bquote("F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
            mtext(text = paste(unique(datum2$species),unique(datum2$individual),sep = "-"),side = 3,outer = T)
            lines(x, y, col="#B22222",lwd=2.5)
            abline(v=round(mean(na.omit(ctmax.half.bs)),2), lty=2,lwd=1.5, col="purple")
            abline(v=round(mean(na.omit(Tcrit.bs)),2), lty=2, lwd=1.5,col="orange")
            points(datum2$fv_fm~datum2$temperature,pch=16,col="#00000095",cex=1.5)
            
            plot(datum2$fv_fm~y3,ylim=c(0,0.88),xlim=c(0,0.88),type="n")
            text(x = 0.03,y = 0.88,labels = bquote(R^2==""))
            text(x = 0.13,y = 0.872,labels = paste(round(summary(model.fit)$r.squared,digits = 3)))
            abline(a = 0,b = 1,col="#808080",lwd=2.5)
            abline(model.fit,col="#B22222",lwd=2.5)
            mtext(text = bquote(Model~"F"[v]*"/"*"F"[m]),side = 1,line = 2,cex = 1)
            mtext(text = bquote(Data~"F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
            points(datum2$fv_fm~y3,pch=16,col="#00000095",cex=1.5)
          }else{# start plotting action 2
            jpeg(filename=paste(graph_dir,#directory
                                as.character(unique(datum2$species)),#species_name
                                as.character(unique(datum2$individual)),#individual_id
                                ".jpg",sep = ""),width = 5500, height = 2500, units = "px", res = 550)
            par(mfrow=c(1,2),oma=c(2,2,2,2), mar=c(2,2,2,2))
            
            plot(datum2$fv_fm~datum2$temperature,ylim=c(0,0.88),xlim=c(24,64),type="n")
            text(x = 55,y = 0.88,labels = bquote("F"[v]*"/"*"F"[m]==""))
            text(x = 61,y = 0.88,labels = paste(round(phia,digits = 3)))
            text(x= 55.8,y = 0.82,labels = bquote("T"[50]==""))
            text(x = 61.1,y = 0.82,labels = paste(round(ctmax[t],digits = 2)))
            text(x= 55.8,y = 0.76,labels = bquote("T"[crit]==""))
            text(x = 61.1,y = 0.76,labels = paste(round(tcrit[t],digits = 2)))
            mtext(text = bquote(Temperature~"("*degree*C*")"),side = 1,line = 2,cex=1)
            mtext(text = bquote("F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
            mtext(text = paste(unique(datum2$species),unique(datum2$individual),sep = "-"),side = 3,outer = T)
            lines(x, y, col="#B22222",lwd=2.5)
            abline(v=round(mean(na.omit(ctmax.half.bs)),2), lty=2,lwd=1.5, col="purple")
            abline(v=round(mean(na.omit(Tcrit.bs)),2), lty=2, lwd=1.5,col="orange")
            points(datum2$fv_fm~datum2$temperature,pch=16,col="#00000095",cex=1.5)
            
            plot(datum2$fv_fm~y3,ylim=c(0,0.88),xlim=c(0,0.88),type="n")
            text(x = 0.03,y = 0.88,labels = bquote(R^2==""))
            text(x = 0.13,y = 0.872,labels = paste(round(summary(model.fit)$r.squared,digits = 3)))
            abline(a = 0,b = 1,col="#808080",lwd=2.5)
            abline(model.fit,col="#B22222",lwd=2.5)
            mtext(text = bquote(Model~"F"[v]*"/"*"F"[m]),side = 1,line = 2,cex = 1)
            mtext(text = bquote(Data~"F"[v]*"/"*"F"[m]),side = 2,line = 2,cex = 1)
            points(datum2$fv_fm~y3,pch=16,col="#00000095",cex=1.5)
            dev.off()
          } # close plotting action 1 & 2
        }# close action to perform if individual has quality data
      }#close individual loop
    
    # 13. Get the output into a dataframe
    TT_output <- data.frame(species = species, ID = ID,# curve control
                            Fv_Fm_Max = Fv_Fm_Max,# maximum fluorescence for an individual
                            ctmax = ctmax,ctmax.uci = ctmax.uci,ctmax.lci = ctmax.lci,# T50
                            tcrit = tcrit,tcrit.uci = tcrit.uci,tcrit.lci = tcrit.lci,# Tcrit
                            R2 = R2)# model fit
    
    Thermal_tol_results[[i]] <- TT_output
    
  }# close species loop
  return(Thermal_tol_results)
}
