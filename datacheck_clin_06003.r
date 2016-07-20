###datacheck.r
##Ch-ch-ch-ch-changes
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()
   
# Set the working directory
   master.dir <- "D:/Hughes/Data"
   scriptname <- "datacheck_clin_06003"
   setwd(master.dir)
   
# Load libraries
   library(ggplot2)
   library(doBy)
   library(Hmisc)
   library(plyr)
   library(grid)
   library(reshape)
   library(stringr)

# Source utility functions file
   source("D:/Hughes/functions_utility.r")
   
# Customize ggplot2 theme - R 2.15.3
   setthemebw2.1()
   
# Organise working and output directories
   working.dir <- paste(master.dir,"RAW_Clinical",sep="/")
   workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

   output.dir <- paste(working.dir,"/",scriptname,"_Output",sep="")
   if(!file.exists(output.dir)){
		dir.create(output.dir)
    }
   
### ------------------------------------- Clinical Data ------------------------------------- ###
### Updated from datacheck_front.r 			#reproducible
   file.name.in <- "RAW_Clinical/rawdata-lena_06003_24hr_allinfo_crcl1.csv"
   file.name.in2 <- "RAW_Clinical/rawdata-lena_06003_24hr_allinfo_CLLn30_average.csv"
   datanew <- read.csv(file.name.in, stringsAsFactors=F, na.strings=c("."))[1:17]
   datasub <- read.csv(file.name.in2,stringsAsFactors=F, na.strings=c("."))[1:16]
   
#------------------------------------------------------------------------------------
#Column names
  #As presented
   names(datanew)
    
  #Sorted
   sort(names(datanew))
  
  #Structure
   str(datanew)

  #datasub the same except missing CrCl1 column
   
#------------------------------------------------------------------------------------
#Plot PK data
  
   with(datanew, table(study, useNA = "always"))
   
   with(datanew, table(dose..mg., useNA = "always"))
   with(datasub, table(dose..mg., useNA = "always"))
   with(datasub, table(dose..mg.,ID))
     
  #Number of patients
  npat <- length(unique(datanew$ID))
  npat
  nsub <- length(unique(datasub$ID))
  nsub

#-------------------------------------------------------------------------------
#Convert datanew to old format 
  
  datanew2 <- data.frame("ID" = datanew$ID, "STUDY" = 6003)
    
  datanew2$STUDY <- 6003
  
  datanew2$XSAMP <- 0
  datanew2$XSAMP[datanew$X.Note=="#"] <- 1
  
  #Group A 	  -> 1 	30 patients	(relapsed CLL)			K. Maddocks et al. 2014
  #Group B 	  -> 2  5 patients (relapsed CLL w/ PR)		K. Maddocks et al. 2014
  #Group C 	  -> 3  24 patients (AML + ALL)				W. Blum et al. 2010
  datanew2$GRP <- 3
  datanew2$GRP[datanew$ID %in% unique(datasub$ID)] <- 1
  datanew2$GRP[datanew$dose..mg.<25&datanew2$GRP!=1] <- 2
  
  #Dose Level 											     		   total patients = 61
  #			Wk1				Wk2				Wk3		
  # ----------------------------------------------------    GROUP A   	npat = 30	
  #	 0	25mg single													   													npat =  3
  #  1	2.5mg daily		5.0mg daily		5.0mg daily					   							npat = 18
  #  2  2.5mg daily		5.0mg daily		7.5mg daily         								npat =  9
  # ----------------------------------------------------    GROUP B	  	npat =  5
  #  1  2.5mg daily   5.0mg daily   5.0mg daily 				   							npat =  3
  #  2  2.5mg daily   5.0mg daily   7.5mg daily 				   							npat =  2  
  # ----------------------------------------------------    GROUP C   	npat = 24
  #  3  25mg daily		25mg daily		25mg daily					  							npat =  4
  #	 4	35mg daily		35mg daily		35mg daily					  							npat =  9
  #  5	50mg daily		50mg daily		50mg daily					  							npat = 10
  #	 6	75mg daily		75mg daily		75mg daily					 								npat =  3
  datanew2$DOSELVL <- 1
  datanew2$DOSELVL[datanew$dose..mg.==25] <- 0
  datanew2$DOSELVL[datanew$dose..mg.==25&datanew2$GRP==3] <- 3
  datanew2$DOSELVL[datanew$dose..mg.==35] <- 4
  datanew2$DOSELVL[datanew$dose..mg.==50] <- 5
  datanew2$DOSELVL[datanew$dose..mg.==75] <- 6
  datanew2$DOSELVL[datanew$ID %in% unique(subset(datanew,dose..mg.==7.5)$ID)] <- 2
  
  datanew2$DOSEMG <- datanew$dose..mg.
        
  datanew2$AMT <- datanew$amt				     #dose in mg
  
  datanew2$RATE <- datanew$rate
   
  datanew2$TIME <- datanew$time..hr.
	
	datanew2$TAD <- datanew$time..hr.
  
  datanew2$WEEK <- ceiling(floor(datanew2$TIME/84)/2)+1
  
  datanew2$DV <- datanew$dv..ug.ml.    			 #ug/ml
 
  datanew2$MDV <- datanew$mdv
    
  #datanew2$BLQ <- 0
  #datanew2$BLQ[datanew$NOTE == "DV_BLQ"] <- 1
  #with(datanew2, table(BLQ, useNA = "always"))
  #datanew2$DV[datanew2$BLQ==1] <- NA
     
  datanew2$LNDV <- log(datanew2$DV)
 
  #datanew2$DNUM <- datanew$Dose.no
  #datanew2$DNUM <- unlist(lapplyBy(~ID, data=datanew2, function(d) impute(d$DNUM)))
  
  #datanew2$OCC <- datanew2$DNUM
  
  datanew2$AGE <- datanew$Age
  
  datanew2$GEND <- datanew$Sex				  	#1 is male, 0 is female
  with(datanew2, table(GEND, useNA = "always"))
  
  datanew2$WT <- datanew$Weight..lbs./2.2		#conversion to kgs
  
  datanew2$HT <- datanew$Height/3.28			#conversion to m	
   
  datanew2$BSA <- 0.007184*datanew2$WT**0.425*datanew2$HT**0.725
  
  datanew2$BMI <- datanew2$WT/datanew2$HT**2
  
    
	#DXCATNUM == 1 -> Chronic Lymphocytic Leukaemia - K. Maddocks et al. 2014
	#DXCATNUM == 2 -> Acute Myeloid Leukaemia		- W. Blum et al. 2010
	#DXCATNUM == 3 -> Acute Lymphoblastic Leukaemia	- W. Blum et al. 2010
  with(datanew, table(Disease, useNA = "always"))
  datanew2$DXCATNUM <- datanew$Disease
  
    #Caucasian 	1
	#??			2
	#??			3
  with(datanew, table(Race, useNA = "always"))
  datanew2$RACE <- datanew$Race
    
  datanew2$RACE2 <- 2
  datanew2$RACE2[datanew$Race == 1] <- 1
  with(datanew2, table(RACE2, useNA = "always"))
  
  datanew2$SECR <- datanew$SeCr..mg.dL.*88.4	#convert from mg/dL to umol/L

 #----------------------------------------------------------------- 
  dataall <- datanew2
   
  dataall <- orderBy(~ID+TIME, data=dataall)
  
#-------------------------------------------------------------------------------
# Check subject numbers
   with(dataall, table(ID))
   
   # Check the dose columns
   with(dataall, table(DOSEMG))
   with(dataall, table(DOSEMG,DOSELVL))	#dose by dose group
 
#----------------------------------------------------------------------------------------------------------------------  
#Calculate dose normalized concentrations

   #Check distribution of DV   
    plotobj <- qplot(x=DV, geom="histogram", data=dataall)
    plotobj
    
    filename.out <- paste(output.dir,"Histogram_DV",sep="/")
    suppressWarnings(to.png(plotobj,filename.out))
    #DV has a wide range of values - 4 orders of magnitude at least
    
   #Check distribution of log DV   
    plotobj <- qplot(x=LNDV, geom="histogram", data=dataall)
    plotobj
    
    filename.out <- paste(output.dir,"Histogram_DVlog",sep="/")
    suppressWarnings(to.png(plotobj,filename.out))
    #DV has a wide range of values - 4 orders of magnitude at least


   #Calculate dose normalised DV
   #Units are ng/ml per mg
  
    dataall$DVNORM <- dataall$DV/dataall$DOSEMG
   
#----------------------------------------------------------------------------------------------------------------------  
#Count missing covariate data
#Missing by Study

   covnames <- as.formula("~AGE+GEND+WT+HT+BSA+BMI+RACE+RACE2+DXCATNUM+SECR")
   covdata <- subset(dataall, select=c("ID","GRP","AGE","GEND","WT","HT","BSA","BMI","RACE","RACE2","DXCATNUM","SECR"))
   
   #Reassign missing
   covdata[covdata==-1] <- NA
  
  #finish off
  missingbystudy <- ddply(covdata, .(GRP), colwise(calculate.percent.missing))
 
  filename.out <- paste(output.dir,"Missing_by_Group.csv",sep="/")
  write.csv(missingbystudy, file=filename.out, row.names=F)
 
  

#Missing by Subject
  missingbysubject <- ddply(covdata, .(GRP,ID), colwise(calculate.percent.missing))
  filename.out <- paste(output.dir,"Missing_by_Subject.csv",sep="/")
  write.csv(missingbysubject, file=filename.out, row.names=F)


#-------------------------------------------------------------------------------
#Subset covariates

   #Keeps missing as -1 - use for categorical summary
   dataallone <- lapplyBy(~ID, data=dataall,  oneperID)
   dataallone <- bind.list(dataallone)
   dim(dataallone)
   
   with(dataallone,table(RATE,useNA="always"))
  
   #Sets missing to NA - use for continuous summary
   covdataone <- lapplyBy(~ID, data=covdata,  oneperID)
   covdataone <- bind.list(covdataone)
   dim(covdataone)
  
  
    dataallone$RACEf <- factor(dataallone$RACE2, labels=c("White or Caucasian","Other"))

    dataallone$SEXf <- as.factor(dataallone$GEND)
    levels(dataallone$SEXf) <- c("female","male")

    dataallone$IDf <- as.factor(dataallone$ID)
  
    dataallone$DXCAT2f <- as.factor(dataallone$DXCATNUM) 
    levels(dataallone$DXCAT2f) <- c("CLL","AML","ALL")

   #DXCAT2 "Relapsed Chronic Lymphocytic Leukaemia" <- 1 
   #DXCAT2 "Acute Myeloid Leukaemia" <- 2 
   #DXCAT2 "Acute Lymphoblastic Leukaemia" <- 3 

#-----------------------------------------------------------------------
#Summary of study characteristics

 #Do all subjects have PK data 
 DVtest <- summaryBy(DV ~ ID, data=dataall, FUN=mean, na.rm=T)
 DVtest <- DVtest[is.na(DVtest$DV.mean)==T,]
 DVtestID <- DVtest$ID
 DVtestID
 #All patients have PK data
 

 #Do all subjects have dose data
  AMTtest <- summaryBy(AMT ~ ID, data=dataall, FUN=mean, na.rm=T)
  AMTtest <- AMTtest[is.na(AMTtest$AMT.mean)==T,]
  AMTtestID <- AMTtest$ID
  AMTtestID
  with(dataallone,table(AMT,useNA="always"))
  #All subjects have dose data
  
 #Do all subjects have nmRATE data
  RATEtest <- summaryBy(RATE ~ ID, data=dataall, FUN=mean, na.rm=T)
  RATEtest <- RATEtest[is.na(RATEtest$RATE.mean)==T,]
  RATEtestID <- RATEtest$ID
  RATEtestID
  with(dataallone,table(RATE,useNA="always"))
  #All subjects have nmRATE data but not in the right place in some cases

 #DV count by Study
 #Calculates data for Report Table 1
 DVcount <- summaryBy(DV ~ GRP, data=dataall, FUN=lengthNA)
 names(DVcount) <- c("Group","DVcount")
 DVcount

 
 #Subject count by Study 
 SUBcount <- ddply(dataall, .(GRP), function(df) count.unique(df$ID))
 names(SUBcount) <- c("Group","SUBcount")
 SUBcount

 #Dose count by Study
 AMTcount <- ddply(dataall, .(GRP), function(df) lengthNA(df$AMT))
 names(AMTcount) <- c("Group","AMTcount")
 AMTcount
 


 #Average DV and AMT per Subject 
  DVsum <- cbind(DVcount,SUBcount,AMTcount)
  DVsum$DVperSUB <-  round(DVsum$DVcount/DVsum$SUBcount,0)
  DVsum$AMTperSUB <-  round(DVsum$AMTcount/DVsum$SUBcount,0)
  DVsum
 
  filename.out <- paste(output.dir,"DVsum.csv",sep="/")
  write.csv(DVsum, file=filename.out)
  



#DV count by Group and DoseLevel
 #Calculates data for Report Table 1
 DVcount <- summaryBy(DV ~ GRP+DOSELVL, data=dataall, FUN=lengthNA)
 names(DVcount) <- c("Group","Dose Level","DVcount")
 DVcount
 
 
#Subject count by Group and DoseLevel
 SUBcount <- ddply(dataall, .(GRP,DOSELVL), function(df) count.unique(df$ID))
 names(SUBcount) <- c("Group","Dose Level","SUBcount")
 SUBcount


#Dose count by Group and DoseLevel
 AMTcount <- ddply(dataall, .(GRP,DOSELVL), function(df) lengthNA(df$AMT))
 names(AMTcount) <- c("Group","Dose Level","AMTcount")
 AMTcount
 

#Average DV and AMT per Subject 
  DVsum <- cbind(DVcount,SUBcount,AMTcount)
  DVsum$DVperSUB <-  round(DVsum$DVcount/DVsum$SUBcount,0)
  DVsum$AMTperSUB <-  round(DVsum$AMTcount/DVsum$SUBcount,0)
  DVsum

  filename.out <- paste(output.dir,"DVsum_DOSELVL.csv",sep="/")
  write.csv(DVsum, file=filename.out)

#DV data present by Group and Week
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0 
  DVcountdata <-  ddply(dataall, .(GRP,ID,WEEK), function(df) DV.present(df$DV))
 

  withDVbyGRPWEEK <- ddply(DVcountdata, .(GRP,WEEK), function(df) sum(df$V1))  #GOLD
  withDVbyGRPWEEK
 

  #Not all subjects with Day 1 data also have Day 4 data
  filename.out <- paste(output.dir,"DVwith_group_week.csv",sep="/")
  write.csv(withDVbyGRPWEEK, file=filename.out)
 
#----------------------------------------------------------------------------------------------------------------------
#Subset some plot data

  plotdata <- subset(dataall)  
  BINnumber <- 3

  plotdata$DOSEMGf <- as.factor(plotdata$DOSEMG)
 
  plotdata$GRPf <- as.factor(plotdata$GRP)
  levels(plotdata$STUDYf) <- paste("Group",levels(plotdata$STUDYf))
  
  plotdata$VOSf <- as.factor(plotdata$VOS)
  levels(plotdata$VOSf) <- c("Placebo","Adrug")

  plotdata$WEEKf <- as.factor(plotdata$WEEK)
  levels(plotdata$WEEKf) <- paste("Week",levels(plotdata$WEEKf))
   
  plotdata$SEXf <- as.factor(plotdata$GEND)
  levels(plotdata$SEXf) <- c("female","male")
  
  plotdata$RACEf <- as.factor(plotdata$RACE2)
  levels(plotdata$RACEf) <- c("White or Caucasian","Other")

  plotdata$DOSELVLf <- as.factor(plotdata$DOSELVL)
  levels(plotdata$DOSELVLf) <- paste("Dose Level",levels(plotdata$DOSELVLf))

  plotdata$DOSE_bin <- cut2(plotdata$DOSEMG, g=BINnumber) 

  plotdata$AGE_bin <- cut2(plotdata$AGE, g=BINnumber) 

  plotdata$WT_bin <- cut2(plotdata$WT, g=BINnumber) 
  
  plotdata$HT_bin <- cut2(plotdata$HT, g=BINnumber)

  plotdata$BSA_bin <- cut2(plotdata$BSA, g=BINnumber) 

  plotdata$BMI_bin <- cut2(plotdata$BMI, g=BINnumber) 

  plotdata$DXCAT2f <- as.factor(plotdata$DXCATNUM)
  levels(plotdata$DXCAT2f) <- c("Relapsed Chronic Lymphocytic Leukaemia","Acute Myeloid Leukaemia","Acute Lymphoblastic Leukaemia")
  
  filename.out <- paste(output.dir,"plotdata.csv",sep="/")
  write.csv(plotdata, file=filename.out, row.names=FALSE)
   
#----------------------------------------------------------------------------------------------------------------------
#Basic PK plots

 #Conc vs Time
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=plotdata)  #, colour=AMTMGf
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=WEEKf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Week")
  plotobj

  filename.out <- paste(output.dir,"Overview_ConcObs_vs_TIME_by_WEEK",sep="/")
  to.png(plotobj,filename.out) 
  
 #Conc vs TIME Week 1
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_DOSELVL",sep="/")
  to.png(plotobj,filename.out) 
  
 #Conc vs TIME Week 1 per ID
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_ID",sep="/")
  to.png(plotobj,filename.out) 
  
   #Conc vs TIME Week 1 per ID (First Half)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(plotdata,ID<36))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_ID_p1",sep="/")
  to.png(plotobj,filename.out) 
  
     #Conc vs TIME Week 1 per ID (Second Half)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(plotdata,ID>35))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_ID_p2",sep="/")
  to.png(plotobj,filename.out) 
  

#----------------------------------------------------------------------------------------------------------------------
#Influence of Covariates

#Function to plot by factor

plotByFactor <- function(factorColname,factorText)
{ 
   
    spanfactor <- 1
  
   #Concentration plots
    plotobj <- NULL
    titletext <- paste("Bdrug Concentrations\n")
    plotobj <- ggplot(data=subset(plotdata,)) 
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y="DV", colour=factorColname), size=2, alpha=0.5)
    #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV"), method=loess, span=spanfactor, se=F, size=1, colour="black")
    plotobj <- plotobj + ggtitle(titletext)                
    plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
    plotobj <- plotobj +  scale_x_continuous("Time after dose (hours)") 
    plotobj <- plotobj + scale_colour_brewer(factorText, palette="Set1")
    plotobj <-  plotobj + facet_wrap(as.formula(paste("~",factorColname,sep="")))  
    plotobj

    filename.out <- paste(output.dir,"/",factorText,"_ConcObs_vs_TAD_facet",sep="")
    to.png.wx2(plotobj,filename.out)

   #Concentration plots
    plotobj <- NULL
    titletext <- paste("Bdrug Concentrations\n")
    plotobj <- ggplot(data=subset(plotdata,))  
    plotobj <- plotobj + geom_point(aes_string(x="TIME", y="DV", colour=factorColname), size=2, alpha=0.5)
    #plotobj <- plotobj + geom_smooth(aes_string(x="TIME", y="DV", colour=factorColname), method=loess, span=spanfactor, se=F, size=1)
    plotobj <- plotobj + ggtitle(titletext)                 
    plotobj <- plotobj + scale_y_log10("Concentration (ng/ml)")
    plotobj <- plotobj + scale_x_continuous("Time after dose (hours)")  
    plotobj <- plotobj + scale_colour_brewer(factorText, palette="Set1")
    plotobj

     filename.out <- paste(output.dir,"/",factorText,"_ConcObs_vs_TAD",sep="")
     to.png.wx2(plotobj,filename.out)
}


#Use the function
  plotByFactor("GRPf","Group")
  plotByFactor("SEXf","SEX")
  plotByFactor("DOSELVLf","Dose Level")
  plotByFactor("WEEKf","Week")
  plotByFactor("DOSE_bin","Binned Dose (mg)")
  plotByFactor("AGE_bin","Binned Age (years)")
  plotByFactor("WT_bin","Binned Weight (kg)")
  plotByFactor("HT_bin","Binned Height (kg)")
  plotByFactor("BSA_bin","Binned BSA (m2)")
  plotByFactor("BMI_bin","Binned BMI (kg per m2)")
  plotByFactor("DXCAT2f","Disease category")


#-------------------------------------------------------------------------
#Summarize Categorical covariates

    dataallone$SUMCOL <- "All Groups"
    
    #Sex
    covCatSexStudy <- with(dataallone,ftable(GRP,SEXf, useNA="ifany", dnn=c("GRP","CATEGORY")))
    covCatSexStudy  <- data.frame("COV"="SEX",covCatSexStudy)


    #Race
    covCatRaceStudy <- with(dataallone,ftable(GRP,RACEf, useNA="ifany", dnn=c("GRP","CATEGORY")))
    covCatRaceStudy  <- data.frame("COV"="RACE",covCatRaceStudy)


    #DXcategory
    covCatDXcatStudy <- with(dataallone,ftable(GRP,DXCAT2f, useNA="ifany", dnn=c("GRP","CATEGORY")))
    covCatDXcatStudy <- data.frame("COV"="DXCAT",covCatDXcatStudy)


    #Collate
    covCatTable <- rbind(covCatSexStudy,covCatRaceStudy,covCatDXcatStudy)
    covCatTable

    #Return to original order
    covCatTable <- orderBy(~COV, covCatTable)  #GOLD - sort by factor levels

    #Define reassignment of column names (if any) rtf allowed
    colNames <- c(
                  "COV","Covariate Code",
                  "CATEGORY","Category",
                   "RACE","Race",
                   "SEX","Gender",
                  "Freq","Count",
				  "DXCAT","Disease"
                  )

    colData <- data.frame(matrix(colNames, byrow=T, ncol=2),stringsAsFactors=F)
    names(colData) <- c("ColNameOld","ColNameNew")

   #Reassign column names
   covCatTable$COV <- gsub.all(colData$ColNameOld,colData$ColNameNew,covCatTable$COV)
   names(covCatTable) <- gsub.all(colData$ColNameOld,colData$ColNameNew,names(covCatTable))

   #Write as an rtf file
   filename.out <- paste(output.dir,"CatCovSummary",sep="/")
   write.csv(covCatTable, file=filename.out)

#-----------------------------------------------------------------------------------------------------
   #Pairs-plot of continuous covariates
   plotobj <- NULL
   covdataonecont <- subset(covdataone, select=c("AGE","WT","BSA","BMI","SECR"))
   plotobj <- ggpairs(na.omit(covdataonecont))
  
      filename.out <- paste(output.dir,"Overview_cont_cov_pairs",sep="/")
      to.png(plotobj,filename.out)


#-----------------------------------------------------------------------------------------------------
 #Pairs-plot of categorical covariates

   covcatdataonecat <- subset(dataallone, select=c("SEXf","RACEf","DXCAT2f"))  

   plotobj <- NULL
   plotobj <- ggpairs(na.omit(covcatdataonecat))

      filename.out <- paste(output.dir,"Overview_cat_cat_pairs",sep="/")
      to.png(plotobj,filename.out)            

      
#-----------------------------------------------------------------------------------------------------------------
#Demographic Summary
#AGE, SEX, WT, BSA

   #Demographics All
   #covDescriptive <- summaryBy(AGE+WT+BSA+BMI~1, data=covdataone, FUN=sumfuncCBIO)
   #covDescriptive <- colwise(formatT)(covDescriptive)
   #filename.out <- paste(output.dir,"demo_summary_all.csv",sep="/")
   #writem.csv(t(covDescriptive), file=filename.out, row.names=F)
    

  #Demographics by GEND
   #covDescriptive <- summaryBy(AGE+WT+BSA+BMI~GEND, data=covdataone, FUN=sumfuncCBIO)
   #covDescriptive <- colwise(formatT)(covDescriptive)
   #filename.out <- paste(output.dir,"demo_summary_GEND.csv",sep="/")
   #writem.csv(t(covDescriptive), file=filename.out, row.names=F)
 

   #Define a custom age bin
   #covdataone$AGEBIN <- cut2(covdataone$AGE, cuts=c(18,50,75,85))
   #covdataone$AGEBIN <- as.numeric(paste(covdataone$AGEBIN))


   #GEND Summary
    #with(covdataone, table(GEND))
  
    #with(covdataone, table(GEND,AGEBIN))
   
 
   #RACE Summary
    #RACEtable <-  with(dataallone, table(RACEf))
        
    #filename.out <- paste(output.dir,"RACEtable.csv",sep="/")
    #write.csv(RACEtable, file=filename.out, row.names=T)
     
    
#-----------------------------------------------------------------------------------------------------------------
#Index plots of covariates

#Customize ggplot2 theme - R 2.15.3
 theme_bw2 <- theme_set(theme_bw(base_size = 20))  
 theme_bw2 <- theme_update(plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
 axis.title.x=element_text(size = 18, vjust = 0),
 axis.title.y=element_text(size = 18, vjust = 1, angle = 90),
 strip.text.x=element_text(size = 14),
 strip.text.y=element_text(size = 14, angle = 90))
 

plotIndexCont <- function(CovColname,CovText)
{ 
  #Debug
  #CovColname <- "HT"
  #CovText <- "Height (cm)"  
  
  plotobj <- ggplot(data=dataallone)
  plotobj <- plotobj + geom_point(aes_string(y=CovColname, x="IDf"), size=3)
  plotobj <- plotobj + scale_x_discrete("Ranked patient index number (ID)")
  CovText <- eval(parse(text=CovText))  #GOLD turn text into expression
  plotobj <- plotobj + scale_y_continuous(name=CovText)
  plotobj <- plotobj + theme(axis.text.x = element_blank()) 
  plotobj

  filename.out <- paste(output.dir,"/IndexPlot_",CovColname,sep="")
  to.png.wx1(plotobj,filename.out)
}


plotIndexCat <- function(CovColname,CovText)
{ 
  #Debug
  #CovColname <- "CYT"
  #CovText <- "Bdrug~Use"  
  
  plotobj <- ggplot(data=dataallone)
  plotobj <- plotobj + geom_point(aes_string(y=CovColname, x="IDf"), size=3)
  plotobj <- plotobj + scale_x_discrete("Ranked patient index number (ID)")
  CovText <- eval(parse(text=CovText))  #GOLD turn text into expression
  plotobj <- plotobj + scale_y_discrete(name=CovText)
  plotobj <- plotobj + theme(axis.text.x = element_blank()) 
  plotobj

  filename.out <- paste(output.dir,"/IndexPlot_",CovColname,sep="")
  to.png.wx1(plotobj,filename.out)
}


#Generate Index plots - CovText is an expression
	plotIndexCont("AGE","Age~(years)")
	plotIndexCont("WT","Weight~(kg)")
	plotIndexCont("BSA","Body~Surface~Area~(m^2)")
	plotIndexCont("BMI","Body~Mass~Index~(kg/m^2)")
	plotIndexCat("GEND","Patient~Sex")
	plotIndexCat("RACE2","Patient~Race")
	plotIndexCat("DXCATNUM","Diagnosis~Category")
	plotIndexCont("SECR","Serum~Creatinine~(umol/L)")

#--------------------
#Data prep
	# [1] "#ID"      "STUDY"    "XSAMP"    "GRP"      "DOSELVL"  "DOSEMG"   "AMT"      "RATE"     "TIME"
	#[10] "TAD"			 "DAY"      "DV"       "MDV"      "LNDV"     "AGE"      "GEND"     "WT"       "HT"
	#[19] "BSA"      "BMI"      "DXCATNUM" "RACE"     "RACE2"    "SECR"     "DVNORM"   "ADDL"     "II"
  dataone <- lapplyBy(~ID, data=dataall,  oneperID)		#one line per point, includes AMT value
  dataone <- bind.list(dataone)
  dataone$RATE <- -2									#fix misplaced rate values
  
  newlines1 <- subset(dataone,DOSELVL==1|DOSELVL==2)
  newlines1$AMT <- 5
  newlines1$DOSEMG <- 5
  newlines1$TIME <- 168
  newlines1$WEEK <- 2
  
  newlines2 <- subset(dataone,DOSELVL==2)
  newlines2$AMT <- 7.5
  newlines2$DOSEMG <- 7.5
  newlines2$TIME <- 336
  newlines2$WEEK <- 3
  
  dataAMT <- rbind(dataone,newlines1,newlines2)
  dataAMT$ADDL <- 6
  dataAMT$ADDL[dataAMT$AMT>20] <- 20
  dataAMT$ADDL[dataAMT$DOSELVL==1&dataAMT$TIME==168] <- 13
  dataAMT$II <- 24
  
  dataDV <- dataall[dataall$TIME>0,]					#all lines not included in dataAMT
  dataDV$AMT <- NA														#remove extra AMT values
  dataDV$RATE <- NA														#fix misplaced rate values
  dataDV$ADDL <- NA										
  dataDV$II <- NA										
  
  dataFIX <- orderBy(~ID+TIME, data=rbind(dataAMT,dataDV))
  colnames(dataFIX)[1] <- "#ID"
  
  colnames(dataFIX)[11] <- "DAY"
	dataFIX$DAY[dataFIX$TIME==0] <- 1
	dataFIX$DAY[dataFIX$TIME!=0] <- ceiling(dataFIX$TIME[dataFIX$TIME!=0]/24)
	dataFIX$TAD <- dataFIX$TIME - (24*(dataFIX$DAY-1))
  
  filename.out <- paste(output.dir,"06003_finaldata.csv",sep="/")
  write.csv(dataFIX, file=filename.out, row.names=FALSE)

