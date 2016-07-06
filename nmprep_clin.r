###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()
   
# Set the working directory
   master.dir <- "D:/Hughes/Data"
   scriptname <- "nmprep_clin"
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
  # ----------------------------------------------------    GROUP A    npat = 30	
  #	 0	25mg single													   npat =  3
  #  1	2.5mg daily		5.0mg daily		5.0mg daily					   npat = 19
  #  2  2.5mg daily		5.0mg daily		7.5mg daily         		   npat =  8
  # ----------------------------------------------------    GROUP B	   npat =  5
  #  1  2.5mg daily     5.0mg daily     5.0mg daily 				   npat =  3
  #  2  2.5mg daily     5.0mg daily     7.5mg daily 				   npat =  2  
  # ----------------------------------------------------    GROUP C    npat = 24
  #  3  25mg daily		25mg daily		25mg daily					   npat =  4
  #	 4	35mg daily		35mg daily		35mg daily					   npat =  9
  #  5	50mg daily		50mg daily		50mg daily					   npat = 10
  #	 6	75mg daily		75mg daily		75mg daily					   npat =  3
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
  
    
	#DXCAT == 1 -> Chronic Lymphocytic Leukaemia - K. Maddocks et al. 2014
	#DXCAT == 2 -> Acute Myeloid Leukaemia		- W. Blum et al. 2010
	#DXCAT == 3 -> Acute Lymphoblastic Leukaemia	- W. Blum et al. 2010
  with(datanew, table(Disease, useNA = "always"))
  datanew2$DXCAT <- datanew$Disease
  
    #Caucasian 	1
	#??			2
	#??			3
  with(datanew, table(Race, useNA = "always"))
  datanew2$RACE <- datanew$Race
    
  datanew2$RACE2 <- 2
  datanew2$RACE2[datanew$Race == 1] <- 1
  with(datanew2, table(RACE2, useNA = "always"))
  
  datanew2$SECR <- datanew$SeCr..mg.dL.*88.4	#convert from mg/dL to umol/L
#-------------------------------------------------------------------------------
#Create nmprep data file

#Week 1 only - Caucasian or Non-Caucasian
#ID TIME AMT EVID DV MDV AGE WT HT GEND RACE SECR DXCAT

  nmcols1 <- datanew2[datanew2$XSAMP==0,-c(2,3,4,5,6,8,10,13,18,19,21)]  #All columns except EVID
  nmcols1$EVID <- 1
  nmcols1$EVID[is.na(nmcols1$AMT)] <- 0

  nmprep1 <- nmcols1[c(1,3,2,13,4,5,6,8,9,7,11,12,10)]
  nmprep1[is.na(nmprep1)] <- "."
  nmprep1
  
  filename.out <- paste(output.dir,"06003LEN_Wk1.csv",sep="/")
  write.csv(nmprep1, file=filename.out)
  
#All Weeks
#ID TIME AMT EVID DV MDV AGE WT HT GEND RACE SECR DXCAT