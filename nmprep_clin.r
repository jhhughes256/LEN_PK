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
  file06003 <- "RAW_Clinical/datacheck_clin_06003_Output/06003_finaldata.csv"
	file05115 <- "RAW_Clinical/datacheck_clin_05115_Output/05115_finaldata.csv"
	file08056 <- "RAW_Clinical/datacheck_clin_08056_Output/08056_finaldata.csv"
	file10016 <- "RAW_Clinical/datacheck_clin_10016_Output/10016_finaldata.csv"
  data06003 <- read.csv(file06003, stringsAsFactors=F)
	data05115 <- read.csv(file05115, stringsAsFactors=F)
	data08056 <- read.csv(file08056, stringsAsFactors=F)
	data10016 <- read.csv(file10016, stringsAsFactors=F)
  
#FINAL datacheck and process for nmprep
  datanew <- rbind(data06003,data05115,data08056,data10016)
	names(datanew)
  str(datanew)
  npat <- length(unique(datanew$X.ID))
  npat
	
	#RATE should be zero for all? Is RATE even necessary as they are all oral doses! 
	with(datanew, table(RATE, useNA = "always"))
	#replace with EVID? far more useful
	colnames(datanew)[8] <- "EVID"
	datanew$EVID <- 0
	datanew$EVID[!is.na(datanew$AMT)] <- 1
	
	#DOSELVL is not uniform between study groups
	print(temp <- with(datanew, table(DOSEMG)))
	#10 different dose levels.. except that patients from 06003 have multiple different dosage levels!
	length(print(names(temp)))
	#There are 9 dose levels +2 unique dose levels from 06003
	with(datanew, table(DOSEMG,STUDY))
	
	dose.levels <- rep(0,length(datanew$DOSELVL))
	dose.levels[datanew$DOSEMG==2.5&datanew$STUDY==8056] <- 1
	dose.levels[datanew$DOSEMG==5&datanew$STUDY==8056] <- 4
	dose.levels[datanew$DOSEMG==15] <- 5
	dose.levels[datanew$DOSEMG==20] <- 6
	dose.levels[datanew$DOSEMG==25] <- 7
	dose.levels[datanew$DOSEMG==30] <- 8
	dose.levels[datanew$DOSEMG==35] <- 9
	dose.levels[datanew$DOSEMG==50] <- 10
	dose.levels[datanew$DOSEMG==75] <- 11
	dose.levels[datanew$GRP<=2&datanew$DOSELVL==1] <- 2
	dose.levels[datanew$GRP<=2&datanew$DOSELVL==2] <- 3
	
	
#Create summary tables
	datanew$DOSELVL <- dose.levels
	
  dataone <- lapplyBy(~X.ID, data=datanew,  oneperID)
  dataone <- bind.list(dataone)
  dim(dataone)
	
	dose.levelsF <- as.factor(dose.levels)
	levels(dose.levelsF) <- c("2.5mg daily","2.5mg daily (escalate to 5.0mg daily)",
														"2.5mg daily (escalate to 7.5mg daily)","5.0mg daily",
														"15mg daily","20mg daily","25mg daily","30mg daily",
														"35mg daily","50mg daily","75mg daily")
	dose.table <- data.frame("Dose Levels" = c(1:11,""),
													 "Regimen" = c(levels(dose.levelsF),"Total Patients"),
													 "No. Patients" = c(with(dataone,table(DOSELVL)),125))
	
	filename.out <- paste(output.dir,"dose_table.csv",sep="/")
  write.csv(dose.table, file=filename.out, row.names=F)

	gt1 <- c(""," ",sort(unique(datanew[datanew$GRP==1,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==2,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==3,]$DOSELVL)),
					 ""," ",sort(unique(datanew[datanew$GRP==4,]$DOSELVL)),
					 ""," ",sort(unique(datanew[datanew$GRP==5,]$DOSELVL)),
					 "  ", #group 6 exists simply to mark lena in comb not a real group
					 ""," ",sort(unique(datanew[datanew$GRP==7,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==8,]$DOSELVL)))
	gt2 <- as.factor(gt1)
	levels(gt2) <- c("STUDY","Group","Group 6 made up of Group 5 patients","2.5mg daily",
									 "50mg daily","75mg daily",
									 "2.5mg daily (escalate to 5.0mg daily)",
									 "2.5mg daily (escalate to 7.5mg daily)",
									 "5.0mg daily","15mg daily","20mg daily",
									 "25mg daily","30mg daily","35mg daily")
	gt3 <- c(with(dataone,table(STUDY))[2],		#specify study patients
					 with(dataone,table(GRP))[1],			#specify group patients
					 temp[1,][temp[1,]!=0],						#specify number of patients for each regimen
					 with(dataone,table(GRP))[2],			#repeat ad nauseam
					 temp[2,][temp[2,]!=0],
					 with(dataone,table(GRP))[3],
					 temp[3,][temp[3,]!=0],
					 with(dataone,table(STUDY))[1],
					 with(dataone,table(GRP))[4],
					 temp[4,][temp[4,]!=0],
					 with(dataone,table(STUDY))[3],
					 with(dataone,table(GRP))[5],
					 temp[5,][temp[5,]!=0],						#group 6 exists simply to mark lena in comb not a real group
					 "(4)",
					 with(dataone,table(STUDY))[4],
					 with(dataone,table(GRP))[6],
					 temp[6,][temp[6,]!=0],						#actually group 7
					 with(dataone,table(GRP))[7],
					 temp[7,][temp[7,]!=0])						#actually group 8
	grp.table <- data.frame("Dose Levels"		= gt1,
													"Regimen" 			= gt2,
													"No. Patients"	= gt3)
													
	#if factor change to character
	i <- sapply(grp.table,is.factor)										
	grp.table[i] <- lapply(grp.table[i],as.character)
	#paste study numbers
  grp.table[which(grp.table$Regimen=="STUDY"),2] <- paste(grp.table[which(grp.table$Regimen=="STUDY"),2],unique(datanew$STUDY))
	#paste group numbers
	grp.table[which(grp.table$Regimen=="Group"),2] <- paste(grp.table[which(grp.table$Regimen=="Group"),2],c(1:5,7,8))
	
	filename.out <- paste(output.dir,"group_table.csv",sep="/")
  write.csv(grp.table, file=filename.out, row.names=F)
	
#Create full nmprep data file	
	filename.out <- paste(output.dir,"fulldata.csv",sep="/")
  write.csv(datanew, file=filename.out, row.names=FALSE)
	
#Prepare nm file
#ID TIME AMT EVID DV MDV ADDL II STUDY GRP DOSELVL AGE GEND WT HT SECR RACE DXCAT

  nmprep <- datanew[c(1,9,7,8,11,12,25,26,2,4,5,14,15,16,17,23,22,20)]
  nmprep[is.na(nmprep)] <- "."
	colnames(nmprep)[c(1,17)] <- c("#ID","RACE")
  
  filename.out <- paste(output.dir,"nmprep_allstudies.csv",sep="/")
  write.csv(nmprep, file=filename.out, quote=FALSE,row.names=FALSE)