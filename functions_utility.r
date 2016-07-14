#--------------------------------------------------------------------------------------------  
# SIMULATION FUNCTIONS
# Simulation function for 2 compartmental model with first order absorption kinetics
   simulate.2comp.abs <- function(ID,AMT,CL,Q,V2,V3,KA,F1)
   {
     k10 <- CL/V2
	 k12 <- Q/V2
	 k21 <- Q/V3
	 apb <- k10+k12+k21            # alpha + beta
	 amb <- k10*k21                # alpha * beta
     alpha <- ((apb)+sqrt((apb)^2-4*amb))/2
	 beta <- ((apb)-sqrt((apb)^2-4*amb))/2
     A <- KA*(k21-alpha)/(V2*(KA-alpha)*(beta-alpha))
	 B <- KA*(k21-beta)/(V2*(KA-beta)*(alpha-beta))
	 Cplasma <- AMT*F1*(A*exp(-alpha*TIME)+B*exp(-beta*TIME)-(A+B)*exp(-KA*TIME))
	 result <- data.frame("TIME"=TIME, "CP"=Cplasma)
	 result
   }   
   
# Process SIM data
  processSIMdata <- function(model.path.file.fit)
   {
   
   # Find simdata and define unwanted lines
      path <- gsub("\\\\", "/", file.choose())
      base.path <- dirname(path)
      #base.path <- dirname(model.path.file.fit)
      setwd(base.path)
      sim.name.in <- basename(path)
      #sim.name.in <- basename(model.path.file.fit)
      sim.name.out <- paste(sim.name.in,".csv", sep="")
      indata <- readLines(sim.name.in)
      tableline <- grep("TABLE NO.  1",indata)  
      headerline <- grep(" ID",indata)          

   # Extract the information in the header line for column names
      header <- indata[headerline[1]]
      header <- scan(textConnection(header), what="char")
      colnum <- length(header)

   # Strip out the unwanted lines
      badlines <- c(tableline,headerline)
      indata <- indata[-badlines]
   
   # Replace spaces for commas and write to a txt file
      indata <- gsub('[[:space:]]+',',',indata)
      writeLines(indata, "SIMtemp.txt")
   
   # Read data as a .csv file and clean it up
      SIMdata <- read.csv("SIMtemp.txt", header=F)
      SIMdata <- SIMdata[-1]
      colnames(SIMdata) <- header
      SIMdata$ID_UNQ <- length(unique(SIMdata$ID))*(SIMdata$SIM-1)+SIMdata$ID
   
   #Write the SIMdata to a file and tidy everything up
      write.csv(SIMdata,sim.name.out, row.names=F)
      file.remove("SIMtemp.txt")
      setwd(master.dir)
   }

#--------------------------------------------------------------------------------------------  
# NON-COMPARTMENTAL FUNCTIONS

# Tmax function (simple)
tmax <- function(x, y)
{ # computes the time of Cmax
 cmax <- max(y)
 tindex <- which(y==cmax)
 tmax <- x[tindex]
 head(tmax, n=1)   #as there can be 2 or more equal Cmax's, choose the first
}

# AUC Function with Linear Trapezoidal Method
 linAUCfunc <- function(dv,time,numobs) 
    {
    
   # Find AUC using trapezoidal method 
   # Area = h*(a+b)/2
   # Area = (t(i+1)-ti)*(Ci+C(i+1))/2
      AUC0t <- sum(diff(time)*(dv[-length(dv)]+dv[-1])/2)
   
      ntail <- 2    # Define ntail to be one value less than desired for first loop as loop adds 1
      bestR2 <- 0   # Define bestR2 so that first loop has something to compare to in first loop
    
	# While number of values being used to determine the slope is less than the number of total values* do the following
	#   - add one to ntail
	#   - define tail values for dv and time
	#   - fit line to log(dv) and time
	#   - define k as the slope and the R2
	#   - if R2 is better than the best R2 so far replace best R2 with new R2 and also replace best k with the new k
	#       *Subtraction of one as C0 may cause an error when logged
      while(ntail<(numobs-1))
      {
        ntail <- ntail+1
        dvtail   <- tail(dv,ntail)
        timetail <- tail(time,ntail)
        fittail  <- lm(log(dvtail) ~ timetail)
   
        k  <- -1*fittail$coefficients["timetail"]
        R2 <- as.numeric(summary(fittail)["r.squared"])
      
        if(R2>bestR2 && k>0)
        {
          bestR2 <- R2
          bestk  <- k
        }
      }
    
      Clast   <- tail(dv,1)     
      AUCtinf <- Clast/bestk     # Calculating extrapolated AUC
      AUC0inf <- AUC0t+AUCtinf   # Adding AUC using trapezoidal method and AUC using extrapolation method
    } 

# AUC Function with Linear Up/Logarithmic Down Trapezoidal Method
   linlogAUCfunc <- function(dv.na,time.na,loq) 
   {
   # Find AUC using trapezoidal method 
   # Define values to be chosen for AUC calculation and give base value for AUC0t
	  n1 <- 1
	  n2 <- 2
	  n3 <- 3
	  n4 <- 4
	  AUC0t <- 0
	  dvtimedf <- na.omit(data.frame(dv.na,time.na))
	  dv <- c(unlist(dvtimedf[1]))
	  time <- c(unlist(dvtimedf[2]))
	  numobs <- length(time)
	  
	  while(n1<numobs)
	  {
	  # Define variables to be used in trapezoidal method
	    c1 <- dv[n1]
		c2 <- dv[n2]
		t1 <- time[n1]
		t2 <- time[n2]
		# Find sum of AUC0-t1 and the new AUCt1-t2
		if(c2>c1)     # if second data point is larger -> Linear Trapezoidal = (t2-t1)*(C1+C2)/2
		{
		  AUCtemp <- (t2-t1)*(c1+c2)/2
		  AUC0t <- sum(AUC0t,AUCtemp)
		}
		if(c2<c1)     # if second data point is smaller -> Logarithmic Trapezoidal = (t2-t1)*(C2-C1)/ln(C2/C1)
		{
		  AUCtemp <- (t2-t1)*(c2-c1)/log(c2/c1)
		  AUC0t <- sum(AUC0t,AUCtemp)
		}
      # Define next values to be chosen for AUC calculation and LOQ searching
		n1 <- n1+1
		n2 <- n2+1
		n3 <- n3+1
		n4 <- n4+1
	  }
    
	# Define R2 to allow comparison, truncate time and dv to remove BLQ values
      ntail <- 3
      bestR2 <- 0	  
      bestk <- 0
	  whichtime <- which(dv==max(dv))              #designate point of cmax
	  adjdv <- dv[(whichtime-1):numobs]            #delete all values before cmax as these are not needed for terminal phase, -1 to allow inclusion of Cmax in regression if required (as this is what WinNonLin does oddly)
	  flag <- 0
  
  	# While number of values being used to determine the slope is less than the number of total values* do the following
    #   - add one to ntail
    #   - define tail values for dv and time
    #   - fit line to log(dv) and time
    #   - define k as the slope and the R2
    #   - if R2 is better than the best R2 so far replace best R2 with new R2 and also replace best k with the new k
    #       *Subtraction of one as C0 may cause an error when logged
	if(length(adjdv)<ntail)
	{ 
	  AUCtinf <- 0
	}else{
      while(ntail<(length(adjdv)) && flag!=1)
      {
		dvtail   <- unlist(tail(adjdv,ntail))
        timetail <- tail(time,ntail)
        fittail  <- lm(log(dvtail) ~ timetail)
   
        k  <- -1*fittail$coefficients["timetail"]
        R2 <- as.numeric(summary(fittail)["r.squared"])
		adjR2 <- 1-((1-R2)*(ntail-1)/(ntail-2))          #adjusted R2 as per WinNonLin user guide
      
        if(adjR2>(bestR2-0.0001) && k>0)                 #if statement as per WinNonLin user guide with added precautions against -ve k vals
        {
          bestR2 <- adjR2
          bestk  <- k
		  ntail <- ntail+1
        }else{
		  flag <- 1
		}
	  }
	  if(bestk==0)
	  {
	    AUCtinf <- 0
	  }else{
	# Calculate AUC from final time to infinite and then add it to AUC from time zero to final time
        Clast   <- tail(dv,1)     
        AUCtinf <- Clast/bestk  
      }		 
	}
    AUC0inf <- AUC0t+AUCtinf  
    } 

# Ke Function
# Could use above to give both AUC and Ke, however this leads to a messy casting and melting situation
kfunc <- function(dv,time,numobs) 
    {
      ntail <- 2    # Define ntail to be one value less than desired for first loop as loop adds 1
      bestR2 <- 0   # Define bestR2 so that first loop has something to compare to in first loop
    
	# While number of values being used to determine the slope is less than the number of total values* do the following
	#   - add one to ntail
	#   - define tail values for dv and time
	#   - fit line to log(dv) and time
	#   - define k as the slope and the R2
	#   - if R2 is better than the best R2 so far replace best R2 with new R2 and also replace best k with the new k
	#       *Subtraction of one as C0 may cause an error when logged
      while(ntail<(numobs-1))
      {
        ntail <- ntail+1
        dvtail   <- tail(dv,ntail)
        timetail <- tail(time,ntail)
        fittail  <- lm(log(dvtail) ~ timetail)
   
        k  <- -1*fittail$coefficients["timetail"]
        R2 <- as.numeric(summary(fittail)["r.squared"])
      
        if(R2>bestR2)
        {
          bestR2 <- R2
          bestk  <- k
        }
      }
      bestk
    } 
	
#--------------------------------------------------------------------------------------------  
# DATA SUMMARY FUNCTIONS

count.unique <- function(x) length(unique(x))

#Convenience function for summarising a dataframe
    headtail <- function(x)
      {
       print(dim(x))
       print(head(x))
       print(tail(x))
      }  

oneperID <- function(x)
{ # returns a single value for each ID of dataframe
  ID <- head(x, n=1)
  ID
}
	  
 #90% confidence interval functions
  CI90lo <- function(x) quantile(x, probs=0.05)
  CI90hi <- function(x) quantile(x, probs=0.95)


#Define a function for length without NA's
   lengthNA <- function(x) length(na.omit(x))
  
 
#Define a function for geometric mean
  geomean <- function(x, na.rm=F)
  {  
  if (na.rm==T) x <- x[is.na(x)==F]
  exp(mean(log(x)))
  }
  #Note x cannot be negative, zero 


#Median and 90% tolerance intervals
  sumfunc90 <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat3 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    stat4 <-  lengthNA(x)
    result <- c("median"=stat1, "low90"=stat2, "hi90"=stat3, "n"=stat4)
    result
  }
  
 #Median and 95% tolerance intervals
  sumfunc95 <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.025, na.rm=T, names=F)  #95%CI
    stat3 <-  quantile(x, probs=0.975, na.rm=T, names=F)
    stat4 <-  lengthNA(x)
    result <- c("median"=stat1, "low95"=stat2, "hi95"=stat3, "n"=stat4)
    result
  }
  
  
#Mean, sd and CV
  sumfuncCV <- function(x)
  {
    stat1 <-  mean(x, na.rm=T)
    stat2 <-  sd(x, na.rm=T)  
    stat3 <-  stat2/stat1*100
    stat4 <-  lengthNA(x)
    result <- c("mean"=stat1, "sd"=stat2, "cv"=stat3, "n"=stat4)
    result
  }
  
  
#Geomean, mean, sd, CV, min, max, n - for Millenium MLN 8237
  sumfuncMLN <- function(x)
  {
    stat1 <-  geomean(x, na.rm=T)
    stat2 <-  mean(x, na.rm=T)
    stat3 <-  sd(x, na.rm=T)  
    stat4 <-  stat3/stat2*100
    stat5 <-  min(x, na.rm=T)  
    stat6 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat7 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    stat8 <-  max(x, na.rm=T) 
    stat9 <-  lengthNA(x)
    result <- c("gmean"=stat1, "mean"=stat2, "sd"=stat3, "cv"=stat4, "min"=stat5, "lo90"=stat6, "hi90"=stat7, "max"=stat8, "n"=stat9)
    result
  }
  
 
#Geomean, mean, sd, CV, min, max, n - for CBIO
  sumfuncCBIO <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  mean(x, na.rm=T)
    stat3 <-  sd(x, na.rm=T)  
    stat4 <-  stat3/stat2*100
    stat5 <-  min(x, na.rm=T)  
    stat6 <-  quantile(x, probs=0.05, na.rm=T, names=F)  #90%CI
    stat7 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    stat8 <-  max(x, na.rm=T) 
    stat9 <-  lengthNA(x)
    result <- c("median"=stat1, "mean"=stat2, "sd"=stat3, "cv"=stat4, "min"=stat5, "lo90"=stat6, "hi90"=stat7, "max"=stat8, "n"=stat9)
    result
  } 
  

#median, mean, sd, CV, 95%CI - for bootstrap parameter summary
  sumfuncBOOT <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  mean(x, na.rm=T)
    stat3 <-  sd(x, na.rm=T)  
    stat4 <-  stat3/stat2*100
    stat5 <-  quantile(x, probs=0.025, na.rm=T, names=F)  #95%CI
    stat6 <-  quantile(x, probs=0.975, na.rm=T, names=F)
    result <- c("median"=stat1, "mean"=stat2, "sd"=stat3, "cv"=stat4, "lo95"=stat5, "hi95"=stat6)
    result
  } 
     
  
 #Summarize distribution by percentiles
  sumfuncPercentile <- function(x)
  {
    stat1 <-  quantile(x, probs=0.05, na.rm=T, names=F) 
    stat2 <-  quantile(x, probs=0.10, na.rm=T, names=F) 
    stat3 <-  quantile(x, probs=0.25, na.rm=T, names=F) 
    stat4 <-  quantile(x, probs=0.50, na.rm=T, names=F) 
    stat5 <-  quantile(x, probs=0.75, na.rm=T, names=F)  
    stat6 <-  quantile(x, probs=0.90, na.rm=T, names=F)  
    stat7 <-  quantile(x, probs=0.95, na.rm=T, names=F)
    result <- c("05perct"=stat1, "10perct"=stat2, "25perct"=stat3, "50perct"=stat4, "75perct"=stat5, "90perct"=stat6, "95perct"=stat7)
    result
  } 

  
#Mean, sd, min and max & n
  sumfuncRange <- function(x)
  {
    stat1 <-  mean(x, na.rm=T)
    stat2 <-  sd(x, na.rm=T)  
    stat3 <-  min(x, na.rm=T)  
    stat4 <-  max(x, na.rm=T)  
    stat5 <-  lengthNA(x)
    result <- c("mean"=stat1,"sd"=stat2, "min"=stat3, "max"=stat4, "n"=stat5)
    result
  }
    
  

#Median etc for boxplot
  sumfuncBOX <- function(x)
  {
    stat1 <-  median(x, na.rm=T)
    stat2 <-  quantile(x, probs=0.025, na.rm=T, names=F) 
    stat3 <-  quantile(x, probs=0.25, na.rm=T, names=F)
    stat4 <-  quantile(x, probs=0.75, na.rm=T, names=F)
    stat5 <-  quantile(x, probs=0.975, na.rm=T, names=F)
    result <- c("median"=stat1, "q025"=stat2, "q25"=stat3, "q75"=stat4, "q975"=stat5)
    result
  }
    


#Define a function for geometric mean and 90% CI of the sem
geomeansemCI <- function(x, na.rm=F)
 #Note x cannot be negative, zero  
  {  
      logx <- log(x)
      logmean <- mean(logx)
      n <- length(x)
      logsem <- sd(logx)/sqrt(n)
        #Critical value of the t-distribution for two one-sided p=0.05
        critt <- qt(.95, df=(n-1)) 
      loglo95 <- logmean - critt*logsem
      loghi95 <- logmean + critt*logsem
      gmean <- exp(logmean)
      glo95 <- exp(loglo95)
      ghi95 <- exp(loghi95)
      result <- c("gmean"=gmean, "glo95"=glo95, "ghi95"=ghi95, "crit.t"=critt)
      result
  }
  
# Percent below the limit of quantification (not including t=0)
   percent.blq <- function(dv,time,blq)
   {
     timedv <- data.frame("TIME"=time,"DV"=dv)
     totaldv <- ifelse(timedv$TIME==0,NA,timedv$DV)
     totaldv <- totaldv[!is.na(totaldv)]
     loqdv <- totaldv[totaldv>=blq]
     percent.bloq <- (1-length(loqdv)/length(totaldv))*100
	 percent.bloq
   }
   
# Tmax
  tmax <- function(dv, time)
{ # computes the time of Cmax
 cmax <- max(dv)
 tindex <- which(dv==cmax)
 tmax <- time[tindex]
 head(tmax, n=1)   #as there can be 2 or more equal Cmax's, choose the first
}

   # Finds 95% CI corrected for mean (this is a ratio, not usable with bioequivalence)
  glohipercent.func <- function(indata,method,variable)
  {
    indata[5] <- indata[3]/indata[2]*100
    indata[6] <- indata[4]/indata[2]*100
    colnames(indata)[5:6] <- c(paste(variable,"_GLO95_PERCENT",sep=""),paste(variable,"_GHI95_PERCENT",sep=""))
    indata <- data.frame("Method"=method,indata)
    indata
   }
   
   # Finds mean and 95% CI for specific dataset
  confint.func <- function(indata,method,variable)
  {
    mean1 <- as.data.frame(sumfunc95(indata[[2]]))
	lo95 <- as.data.frame(sumfunc95(indata[[3]]))
    hi95 <- as.data.frame(sumfunc95(indata[[4]]))
    VARIABLE <- c(paste(variable,"_LO95",sep=""),paste(variable,"_MEAN",sep=""),paste(variable,"_HI95",sep=""))
	ANALYSIS <- method
	CILO95 <- c(lo95[2,], mean1[2,], hi95[2,])
	MEAN <- c(lo95[1,], mean1[1,], hi95[1,])
	CIHI95 <- c(lo95[3,], mean1[3,], hi95[3,])
	all <- data.frame(VARIABLE,ANALYSIS,CILO95,MEAN,CIHI95)
	all
  }
  
  # Determine bioequivalence
  bioq.func <- function(outputdf,limitlo,limithi,ctl.name)
  {
    outputdf$PF_FREL <- ifelse(outputdf$FREL_GLO95 < limitlo | outputdf$FREL_GHI95 > limithi,1,0)
    probtable <- ddply(outputdf, .(SIM_ID), function(df) CalcProb(df$PF_FREL))
    probtable <- data.frame("Metric"="FREL","Data"=ctl.name, probtable)
  }
  
    crat.func <- function(outputdf,limitlo,limithi,ctl.name)
  {
    outputdf$PF_CMAX <- ifelse(outputdf$CMAX_GLO95 < limitlo | outputdf$CMAX_GHI95 > limithi,1,0)
    probtable <- ddply(outputdf, .(SIM_ID), function(df) CalcProb(df$PF_CMAX))
    probtable <- data.frame("Metric"="CRAT","Data"=ctl.name, probtable)
  }
  
     percent.blq.time <- function(df)
   {
     subtimes <- as.numeric(unlist(unique(df[1])))
	 ntimes <- length(subtimes)
	 percent.blq <- 0
	 for(i in 1:ntimes)
	 {
	   T <- subset(df,df[1]==subtimes[i])
	   percent.blq <- c(percent.blq,mean(as.numeric(unlist(T[4]))))
	 }
	 percent.blq <- percent.blq[-1]
	 percent.blq
   }
  
  
   # Works for specific ggplot, finds mean of a number for binned times
   percent.blq.time.multi <- function(df)
   {
     subtimes <- as.numeric(unlist(unique(df[1])))
	 ntimes <- length(subtimes)
	 percent.blq1 <- 0
	 percent.blq2 <- 0
	 form <- 0
	 for(i in 1:ntimes)
	 {
	   T <- subset(df,df[1]==subtimes[i])
	   T1 <- subset(T,T[2]=="Formulation 1")
	   T2 <- subset(T,T[2]=="Formulation 2")
	   percent.blq1 <- c(percent.blq1,mean(as.numeric(unlist(T1[4]))))
	   percent.blq2 <- c(percent.blq2,mean(as.numeric(unlist(T2[4]))))
	 }
	 percent.blq1 <- percent.blq1[-1]
	 percent.blq2 <- percent.blq2[-1]
	 percent.blq <- c(percent.blq1,percent.blq2)
	 form <- rep(c("Formulation 1","Formulation 2"),each=14)
	 tab <- data.frame(form,percent.blq)
	 tab
   }
   
  ### Assess Bioequivalence
 # Assign Pass/Fail flag to Confidence Intervals
   # 0 is CI within limits, 1 is CI outside limits
   CalcProb <- function(x)
   # Probability of 0 for binary events coded as 0 and 1
    {
     prob <- sum(x)/length(x)
     prob <- 1-prob
     c("p"=prob)  #p of zero
    }
	
# Checking for type1 and type2 error against reference data
  # 1 is positive, 0 is negative
	  errortype.func <- function(ref,test)   # type=1 -> Type1 Error  type=2 -> Type2 Error
   {
	 T1error <- gsub(TRUE,1,ref>test)	 # REF>TEST <- Type1 Error
	 T1error <- gsub(FALSE,0,T1error)
	 T1error <- gsub("NA",0,T1error)
	 T2error <- gsub(TRUE,1,ref<test)    # REF<TEST <- Type2 Error
	 T2error <- gsub(FALSE,0,T2error)
	 T2error <- gsub("NA",0,T2error)
	 data.frame("pT1"=T1error,"pT2"=T2error)
   }
      # Finds true SIMID for bioqtables, of particular use in M3 where NONMEM can fail to create fit file
	  # Then uses errortype.func to find type1 and type2 error
      errortype.process <- function(test,ref,fitfail,nsim,tag)
   {
     nfail <- length(fitfail)
	 trueSIMID <- (1:nsim)[-c(fitfail)]
     if(fitfail[1]!=0)
     {
       test$SIM_ID <- trueSIMID
	   missingdf <- data.frame(Metric=rep("FREL",times=nfail),"Data"=rep(tag,times=nfail),"SIM_ID"=fitfail,"p"=rep("NA",times=nfail))
	   truetable <- rbind(test,missingdf)
	   truetable <- orderBy(~SIM_ID,truetable)
	   result <- truetable
     }else{
	   result <- test
	 }
	 final <- data.frame(result,errortype.func(ref$p,result$p))
	 final
   }

### General R Functions	
# Allows combination of two or more vectors such that the first occurrence of a non NA is taken
# > a <- c(1,  2,  NA, NA, NA)
# > b <- c(NA, NA, NA, NA, 6)
# > c <- c(7,  8,  NA, 9, 10)
# > coalesce2(a,b,c)
# # [1]  1  2 NA  9  6
   coalesce2 <- function(...) {
      Reduce(function(x, y) {
         i <- which(is.na(x))
         x[i] <- y[i]
         x},
      list(...))
   }

# Control Stream creating function 
# Can be used for mass replacement of specific parts of a .txt file
# pat and repl are vectors that represent each pattern to be replaced and each replacement respectively
# ctl.rep
   
   gsub.all <- function(pat,repl,x)
   {
     vec <- x
     if(length(pat)==length(repl))
	 {
       for(i in 1:length(pat))
	   {
	     vec <- gsub(pat[i],repl[i],vec,fixed=TRUE)
	   }
	   vec
	 }else{
	   warning("length(pattern)!=length(replacement): Amount of values to be replaced is not equal to amount of values given")
	 }
   }
   
   #Utility function to bind a list of dataframes - must have matching columns
  bind.list <- function(x)
   {
    #Bind a list of smaller dataframes into one big dataframe
    #Access the 2nd level of the list with [[x]]
   alldata <- x[[1]]
   for (i in 2:length(x))
    {
    alldata <- rbind( alldata, x[[i]] )
    as.data.frame(alldata)
    }
   alldata 
   }
   
   
   
# NONMEM Batch File creating function (also directs NONMEM to the .csv file by changing "dataname" in the controls stream)

   nm.prep <- function(nsim,ncore,nid,EST.file,ctl.name,SIM.name.out,limdata,limobs,blq,ctlref,nmbat,amt)
   {
   for (i in 1:nsim)
   {
	 file.name <- paste(EST.file,ctl.name,"_model",i,".ctl",sep="")
	 data.name <- paste(SIM.name.out,"_model",i,".csv",sep="")
     tempctl <- sub("dataname",data.name,ctlref)
	 writeLines(tempctl,file.name)
# Create .csv for NONMEM input for each SIM
     subdata <- subset(limdata,SIM==i)
	 SIM <- rep(i,times=2*limobs*nid)
     EVID <- ifelse(subdata$TIME==0,4,0) 
     BLQ <- ifelse(subdata$DV<blq,1,0)
	 BLQ <- ifelse(subdata$TIME==0,0,BLQ)  
     AMT <- ifelse(EVID==4,amt,".")
	 MDV <- ifelse(AMT==amt,1,0)
	 MDV <- ifelse(subdata$DV<blq,1,MDV)
     nlmedata <- data.frame(subdata$ID,subdata$TIME,AMT,EVID,subdata$FORM,subdata$DV,MDV,SIM,BLQ)
     colnames(nlmedata) <- c("#ID","TIME","AMT","EVID","FORM","DV","MDVX","SIM","BLQ")
     file.name <- paste(EST.file,"_model",i,".csv",sep="")
     write.csv(nlmedata,file.name,quote=FALSE,row.names=FALSE)
# Create command lines to be split into seperate .bats
     tempbat <- paste("call nmgo ",SIM.name.out,ctl.name,"_model",i,".ctl",sep="")
	 nmbat[i] <- tempbat
   }
   nmbat
   }
   
   # Gives thetas found by NONMEM
   nlme.fit <- function(SIM.name.out,FIT.dir,EST.dir,ctl.name,nsim,nid,limobs)
   {
     fit2sim <- data.frame(matrix(NA, nrow = nsub*2, ncol = 11))
     colnames(fit2sim) <- c("UNQ_ID","STUD_ID","SIM_ID","AMT","F1","CL","Q","V2","V3","KA","AUC")
	 ctl.num <- ifelse(ctl.name=="M1",1,3)
	 rnum <- 1
	 fail.fit <- 0
     for (i in 1:nsim)
	 {
       fit.name.in <- paste(SIM.name.out,ctl.num,"_model",i,sep="")
	   fit.dir.out <- paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep="")
	   fit.file <- paste(EST.dir,"/",fit.name.in,".nm7/",fit.name.in,".fit",sep="")
	   if(file.exists(fit.file))                                                     # Due to M3 failing to create a .fit file in ~0.5% of cases
	   {
	     fitdata <- read.table(file=fit.file, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
	     fitdata <- cbind(rep(1:nid+nid*(i-1),each=limobs*2),fitdata)
	     colnames(fitdata)[1] <- "UNQ_ID"
	     write.csv(fitdata, file=fit.dir.out, row.names=FALSE) 
	     fit.temp <- read.csv(paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep=""))
	     fit.temp <- subset(fit.temp,fit.temp$TIME==0)
	     fit.temp <- fit.temp[-c(4,6,7,8,9,17,18,19,20,21)]
		 fit.temp[3] <- rep(rnum,times=nid*2)
	     fit2sim[(1:(nid*2))+nid*2*(rnum-1),] <- fit.temp
		 rnum <- rnum+1                                                              # Record how many sims were successful to ensure safe use of data with universal functions
	   }else{
	     fail.fit <- c(fail.fit,i)                                                   # Record rows that did not give a .fit file
	   }
     }
	 if(length(fail.fit)>1)
	 {
	   fail.fit <- fail.fit[-1]                                                      # If no sims failed to give .fit file, output==0, if sims did fail to give .fit file, remove 0 from output
	 }
     write.csv(na.omit(fit2sim), file=paste(SIM.file,ctl.name,"NMTHETAS.csv", sep="_"), row.names=FALSE)
	 c(rnum-1,fail.fit)
   }
   
   nlme.fit2 <- function(SIM.name.out,FIT.dir,EST.dir,ctl.name,nsim,nid,limobs) #same as above but doesnt change files
   {
     fit2sim <- data.frame(matrix(NA, nrow = nsub*2, ncol = 11))
     colnames(fit2sim) <- c("UNQ_ID","STUD_ID","SIM_ID","AMT","F1","CL","Q","V2","V3","KA","AUC")
	 ctl.num <- ifelse(ctl.name=="M1",1,3)
	 rnum <- 1
	 fail.fit <- 0
     for (i in 1:nsim)
	 {
       fit.name.in <- paste(SIM.name.out,ctl.num,"_model",i,sep="")
	   fit.dir.out <- paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep="")
	   fit.file <- paste(EST.dir,"/",fit.name.in,".nm7/",fit.name.in,".fit",sep="")
	   if(file.exists(fit.file))                                                     # Due to M3 failing to create a .fit file in ~0.5% of cases
	   {
	     fitdata <- read.table(file=fit.file, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
	     fitdata <- cbind(rep(1:nid+nid*(i-1),each=limobs*2),fitdata)
	     colnames(fitdata)[1] <- "UNQ_ID"
	     write.csv(fitdata, file=fit.dir.out, row.names=FALSE) 
	     fit.temp <- read.csv(paste(FIT.dir,"/",ctl.name,"_fit",rnum,".csv",sep=""))
	     fit.temp <- subset(fit.temp,fit.temp$TIME==0)
	     fit.temp <- fit.temp[-c(4,6,7,8,9,17,18,19,20,21)]
		 fit.temp[3] <- rep(rnum,times=nid*2)
	     fit2sim[(1:(nid*2))+nid*2*(rnum-1),] <- fit.temp
		 rnum <- rnum+1                                                              # Record how many sims were successful to ensure safe use of data with universal functions
	   }else{
	     fail.fit <- c(fail.fit,i)                                                   # Record rows that did not give a .fit file
	   }
     }
	 if(length(fail.fit)>1)
	 {
	   fail.fit <- fail.fit[-1]                                                      # If no sims failed to give .fit file, output==0, if sims did fail to give .fit file, remove 0 from output
	 }
	 c(rnum-1,fail.fit)
   }
 
# GGPLOT FUNCTIONS 
#--------------------------------------------------------------------------------------------
# Set the theme
   setthemebw2 <- function()
   {
      theme_bw2 <- theme_set(theme_bw(base_size = 20))  
      theme_bw2 <- theme_update(plot.margin = unit(c(1.1,1.1,3,1.1), "lines"),
      axis.title.x=element_text(size = 18, vjust = 0),
      axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
      strip.text.x=element_text(size = 14),
      strip.text.y=element_text(size = 14, angle = 90))
   }
   setthemebw2.1 <- function()
   {
    theme_bw2 <- theme_set(theme_bw(base_size = 20))  
    theme_bw2 <- theme_update(plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
    axis.title.x=element_text(size = 18, vjust = 0),
    axis.title.y=element_text(size = 18, vjust = 1, angle = 90),
    strip.text.x=element_text(size = 14),
    strip.text.y=element_text(size = 14, angle = 90))
   }
   
#--------------------------------------------------------------------------------------------  
# to.png functions

to.png <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=900, height=700, pointsize=pointsizein)
  print(plotobj.in)
  
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8))
		 
  dev.off()
 }
 
to.png.sqr <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=650, height=500, pointsize=pointsizein)
   print(plotobj.in)
   
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
  dev.off()
 }    


to.png.tiny <- function(plotobj.in,maintext.in,pointsizein=12)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=250, height=250, pointsize=pointsizein)
   print(plotobj.in)
   
  # #Stamp the plot
  # time.date <- Sys.time()
  # working.dir <- getwd()
  # grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            # x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            # just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
   dev.off()
 }      

 
 to.png.wx1 <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=700, height=700, pointsize=pointsizein)
   print(plotobj.in)
   
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
  dev.off()
 }      
 
 to.png.wx2 <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=1350, height=700, pointsize=pointsizein)
   print(plotobj.in)
   
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8)) 
   
  dev.off()
 }      

to.png.hx1.5 <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=900, height=975, pointsize=pointsizein)
  print(plotobj.in)
  
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8))
		 
  dev.off()
 }  
 
to.png.hx2 <- function(plotobj.in,maintext.in,pointsizein=14)
 {
  png.file.name <- paste(maintext.in,".png",sep="")
  png(file=png.file.name, width=900, height=1200, pointsize=pointsizein)
  print(plotobj.in)
  
  #Stamp the plot
  time.date <- Sys.time()
  working.dir <- getwd()
  grid.text(paste("Source: ",working.dir,maintext.in,time.date),
            x = unit(0.025,"npc"), y = unit(0.015, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 8))
		 
  dev.off()
 } 
 
 runaov2 <- function(df, USE)
 {
   #debug
   #df <- EXP.data[EXP.data$REP==1,]
   #USE <- "AUC"
   #uselog <- T
   
   #BE aov for study replicate i
   df$metric <- df[,USE]
   result <- aov(log(metric)~FORM, data=df) 
   result2 <- summary(result)
   
   #This function avoids changing the contrasts for the aov 
   aovtable <- model.tables(result,"means", se=T)
   int <- confint(result, level=0.9) # FUNCTION TO CALCULATE CONFIDENCE INTERVAL OF LN DATA
   
   #Extracting aov results
   Xt <- aovtable$tables$FORM[2]
   Xr <- aovtable$tables$FORM[1]
   
   pointestimate <- exp(Xt-Xr)
   pointestimate
   
   CI90lovalue <- exp(int[2,1])
   CI90lovalue
   
   CI90hivalue <- exp(int[2,2])
   CI90hivalue
   
   #Bioequivalence test 1=bioequivalent
   bioflag <- 0
   if (CI90lovalue > 0.8 & CI90hivalue < 1.25) bioflag <- 1
   
   BEresultsi <- data.frame("Metric"=USE,"pointestimate"=pointestimate,"lowerCI"=CI90lovalue,"upperCI"=CI90hivalue,"BE"=bioflag)
   
 }
 
       errortype.process2 <- function(test,ref,fitfail,nsim,tag)
   {
     test2 <- data.frame(test[,1],test[,2],test[,6])
	 colnames(test2) <- c("SIM","Metric","BE")
     if(fitfail[1]!=0)
     {
	   nfail <- length(fitfail)
	   trueSIMID <- (1:nsim)[-c(fitfail)]
       test2[1] <- trueSIMID
	   missingdf <- data.frame("SIM"=fitfail,"Metric"=rep("AUC",times=nfail),"BE"=rep("NA",times=nfail))
	   truetable <- rbind(test2,missingdf)
	   truetable <- orderBy(~SIM,truetable)
	   result <- truetable
     }else{
	   result <- test2
	 }
	 final <- errortype.func2(ref$BE,result$BE)
	 final
   }
   
   	  errortype.func2 <- function(ref,test)   
   {
	 T1error <- gsub(TRUE,1,ref<test)	 
	 T1error <- gsub(FALSE,0,T1error)
	 T1error <- gsub("NA",0,T1error)
	 T2error <- gsub(TRUE,1,ref>test)
	 T2error <- gsub(FALSE,0,T2error)
	 T2error <- gsub("NA",0,T2error)
	 data.frame("pT1"=T1error,"pT2"=T2error)
   }
   
   ################################################################################
# Some functions used via R2HTML to comment an R script and produce html output
# needs library(R2HTML)
# Very simple, uncluttered way of converting R script into commented HTML

#HTML command functions are:
#HE(n) - echos the last "n" R commands and writes them to an HTML file
#HEO(n) - echos the last "n" R commands and their output and writes them to an HTML file
#HC(Hlevel,text) - writes a series of comments to an HTML file - heading level is specified by Hlevel, the text is a string
#HG(file.name) - inserts a jpeg graph file into the html file
#HGG2(plotobj,maintext) - creates a jpeg file of a ggplot2 object, then inserts the jpeg into the html file - maintext becomes the both the name of the graph and title of the file
#HT(file.name) - retrieves a *.csv file as a dataframe, then writes the dataframe as an html table
#HL(file.name) - inserts a hyperlink to a file into the html file


#Initiate the R2HTML file
#CSS.file <- ("R2HTML.css")
#HTMLInitFile(outdir=working.dir, filename="example",  HTMLframe=F, CSSFile=CSS.file)  #don't add extension

#code goes here
#HC(1,"PAGANZ 10 Workshop - Predictive Checks using in NM7, WfN and R")
#HC(2,"Richard Upton")
#HC(3,"This example illustrates the power of using a scripting language like R")
#HC(5,"The analysis draws on the principles of reproducible research, as discussed in this paper by Gentleman and Lang:")
#HL("Gentleman_2004_Reproducible_Research.pdf")

#Some html tags are allowed, e.g <br>  <b>Bold</b>

#HTMLEndFile(file="example.html")   #add extension

#Standard R2HTML commands can also be used:
 #exampledata <- read.csv("testdata_ALBWTcov2.csv")
 #HTML(exampledata)


#Echos and command only in the HTML file!
HE <- function(num.commands)
#HTML ECHO of command - echos the last "n" R commands and writes them to an HTML file
  {
   savehistory()
   command.history <- readLines(".Rhistory")
   n <- num.commands
   l <- length(command.history)
   last.commands <- command.history[(l-n):(l)]
   for (i in 1:n)
   {
   HTML(paste(">",last.commands[i],sep=""))
   HTML("")
   }
  }

#Echos and command AND output in the HTML file - brilliant!
HEO <- function(num.commands)
#HTML ECHO of command - echos the last "n" R commands and writes them to an HTML file
  {
   savehistory()
   command.history <- readLines(".Rhistory")
   n <- num.commands
   l <- length(command.history)
   last.commands <- command.history[(l-n):(l)]
   for (i in 1:n)
   {
   HTML(paste(">",last.commands[i],sep=""))
   HTML((eval(parse(text=last.commands[i]))))
   HTML("")
   }
  }

HC <- function(Hlevel,text)
#HTML COMMENT - writes a series of comments to an HTML file
#The HTML heading level is specified by Hlevel, the text is a string
#html code can be inserted and will (probably!) be understood
  {
  if (Hlevel == 1)
  HTML(paste("<H1>",text,"</H1>"))
  if (Hlevel == 2)
  HTML(paste("<H2>",text,"</H2>"))
  if (Hlevel == 3)
  HTML(paste("<H3>",text,"</H3>"))
  if (Hlevel == 4)
  HTML(paste("<H4>",text,"</H4>"))
  if (Hlevel == 5)
  HTML(paste("<p>",text,"</p>"))
  }


HG <- function(file.path)
#HTML GRAPH - inserts a png graph file into the html file
#requires master.dir from wfnviaR settings file
 {
  #debug file.path <- "1_1comp_ka_lag.g77/conc_vs_time.jpg"
  HTMLInsertGraph(GraphFileName=file.path, WidthHTML=600)
 }


HGG2 <- function(plotobj,maintext)
#HTML GGPLOT2 GRAPH - creates a png file of a ggplot2 object, then inserts the png into the html file
#Optimised for metadata style script - see P:\MLN8237PKupdate\Datacheck6\datacheck6.html_script.R
 {
 png.file.name <- paste(maintext,".png",sep="")
 to.png(plotobj,maintext)
 HG(paste("file://",png.file.name,sep=""))
 }
 

HGG2long <- function(plotobj,maintext)
#HTML GGPLOT2 GRAPH - creates a png file of a ggplot2 object, then inserts the png into the html file
#Optimised for metadata style script - see P:\MLN8237PKupdate\Datacheck6\datacheck6.html_script.R
 {
 png.file.name <- paste(maintext,".png",sep="")
 to.png.long(plotobj,maintext)
 HG(paste("file://",png.file.name,sep=""))
 }


HGG2wide <- function(plotobj,maintext)
#HTML GGPLOT2 GRAPH - creates a png file of a ggplot2 object, then inserts the png into the html file
#Optimised for metadata style script - see P:\MLN8237PKupdate\Datacheck6\datacheck6.html_script.R
 {
 png.file.name <- paste(maintext,".png",sep="")
 to.png.wide(plotobj,maintext)
 HG(paste("file://",png.file.name,sep=""))
 }
 

HT <- function(file.path)
#HTML TABLE - retrieves a *.csv file as a dataframe, then writes the dataframe as an html table
 {
  #debug file.path <- "1_1comp_ka_lag.g77/1_1comp_ka_lag.fit.param.csv"
  temp.data.frame <- read.csv(file.path, stringsAsFactors=F)
  print(temp.data.frame)                                      #output in console
  HTML(temp.data.frame)                                #output in html
  rm(temp.data.frame)
 }



HL <- function(file.path)
#HTML LINK - inserts a hyperlink to a file into the html file
#requires master.dir from wfnviaR settings file
 {
 #Debug file.path <- "concdata.csv"
 #link.file <- paste(master.dir,file.path, sep="/")
 #link.file <- paste(master.dir,file.path, sep="/")
 #paste("<a href=",file.path,">",link.file,"</a>", sep=" ")
 link.string <- paste("<a href=",file.path,">",file.path,"</a>", sep=" ")
 HTML(link.string)
 #<a href="../rawdata1.xls">../rawdata1.xls</a>
 
 #HL("concdata.csv")
 }
 
 
  #------------------------------------------------------------------------------------
#Function for filling in gaps in covariate data - for one ID - use doBy to apply to all ID
#Only applies for covariates that don't change with time
   #test0 <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05)
   #test1 <- c(NA,0.05,0.05,NA,0.05,NA,0.05)
   #test2 <- c(NA,"2_0.0833_0.05_1","2_0.0833_0.05_1",NA,"2_0.0833_0.05_1",NA,"2_0.0833_0.05_1")
   #test3 <- c(NA,NA,NA,NA)
   #test4 <- c(0.05,0.01,0.05,0.05,NA,0.05,0.05)
   
   fill.missing <- function(x)   #GOLD
    {
    all.fill <- unique(x[is.na(x)==F])  #All non-missing values
    fill <- all.fill[1]  #value to replace missing values
       if (length(all.fill)>1) stop("Can't fill with non-unique values") else x[is.na(x)==T] <- fill
    x
    }
    
    #fill.missing(test0)
    #fill.missing(test1)
    #fill.missing(test2)
    #fill.missing(test3)
    #fill.missing(test4)

#----------------------------------------------------------------------------------------------  

###Imputation functions###

locf <- function (x)
#Last observation carried forward
#Finds an NA and carries forward the previous value  
  {
    good <- !is.na(x)
    positions <- seq(length(x))
    good.positions <- good * positions
    last.good.position <- cummax(good.positions)
    last.good.position[last.good.position == 0] <- NA
    x[last.good.position]
}                                      


nocb <- function (x)
#Next observation carried backward
#Reverses the locf function  
  {
   rev(locf(rev(x)))
  }  

 
impute <- function(x)
#Function that runs locf first, then nocb
{
    x <- locf(x)
    nocb(x)
    x
}    

#-------------------------------------------------------------------------------------------
###Functions for counting missing data###

#see also lengthNA

#Function for calculating percent missing - use apply to do column by column  
 calculate.percent.missing <- function(x)
    {
    length(x[is.na(x)==T])/length(x)*100
    } 

#Function returns T if any NA's
 anyNA <- function(x) any(is.na(x))

#Function returns T if all NA's
 allNA <- function(x) all(is.na(x))

# #Testing
# Concentration <- c(1,2,3,4)
# percent.missing(Concentration)
# anyNA(Concentration)
# allNA(Concentration)

# Concentration <- c(1,2,3,NA,4)
# percent.missing(Concentration)
# anyNA(Concentration)
# allNA(Concentration)

# Concentration <- c(NA,NA,NA,NA,NA)
# percent.missing(Concentration)
# anyNA(Concentration)
# allNA(Concentration)

# #Run this at STUDY level and ID level

#-----------------------------------------------------------------------------
#Functions for listing and finding objects

# #Search for objects in workspace with wildcards
# ls(pattern=glob2rx("final*"))

# #Show all objects
# ls.str()
# print(ls.str(), max.level = 0)# don't show details


# #Functions
# #how do the functions look like which I am using?
# lsf.str() 

# #write the function to a text file
# writeLines(lsf.str(), "function_list.txt")

# #list function by wildcard
# lsf.str(pattern = "Draw")


# #Dataframes
# #dataframes are lists in this context
# ls.str(mode = "list") 

# print(ls.str(mode = "list"), max.level=0) 

# print(ls.str(mode = "list",pattern=glob2rx("*data*")), max.level=0)

# writeLines(print(ls.str(mode = "list",pattern=glob2rx("*data*")), max.level=0), "dataframe_list.txt")


# #show packages
# search() 


#------------------------------------------------------------------------------------------------------ 
#Alternatively - save all functions in an environment to text
#from the interweb!

save.functions.from.env <- function(file = "d:\\temp.r")
{
    # This function will go through all your defined functions and write them into "d:\\temp.r"
    # let's get all the functions from the envoirnement:
    funs <- Filter(is.function, sapply(ls( ".GlobalEnv"), get))
 
    # Let's 
    for(i in seq_along(funs))
    {
        cat(    # number the function we are about to add
            paste("\n" , "#------ Function number ", i , "-----------------------------------" ,"\n"),
            append = T, file = file
            )
 
        cat(    # print the function into the file
            paste(names(funs)[i] , "<-", paste(capture.output(funs[[i]]), collapse = "\n"), collapse = "\n"),
            append = T, file = file
            )
 
        cat(
            paste("\n" , "#-----------------------------------------" ,"\n"),
            append = T, file = file
            )
    }
 
    cat( # writing at the end of the file how many new functions where added to it
        paste("# A total of ", length(funs), " Functions where written into", file),
        append = T, file = file
        )
    print(paste("A total of ", length(funs), " Functions where written into", file))
}
 
 #save.functions.from.env(file="text.txt") # this is how you run it

repeat.before = function(x) {   # repeats the last non NA value. Keeps leading NA
    ind = which(!is.na(x))      # get positions of nonmissing values
    if(is.na(x[1]))             # if it begins with a missing, add the 
          ind = c(1,ind)        # first position to the indices
    rep(x[ind], times = diff(   # repeat the values at these indices
       c(ind, length(x) + 1) )) # diffing the indices + length yields how often 
}                               # they need to be repeated
 