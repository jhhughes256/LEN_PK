###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
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
  source("E:/Hughes/functions_utility.r")

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
	dose.levels[datanew$DOSEMG==7.5&datanew$STUDY==8056] <- 5
	dose.levels[datanew$DOSEMG==15] <- 6
	dose.levels[datanew$DOSEMG==20] <- 7
	dose.levels[datanew$DOSEMG==25] <- 8
	dose.levels[datanew$DOSEMG==30] <- 9
	dose.levels[datanew$DOSEMG==35] <- 10
	dose.levels[datanew$DOSEMG==50] <- 11
	dose.levels[datanew$DOSEMG==75] <- 12
	dose.levels[datanew$GRP<=2&datanew$DOSELVL==1] <- 2
	dose.levels[datanew$GRP<=2&datanew$DOSELVL==2] <- 3

	temp <- with(datanew,table(MDV,X.ID))
	#ID with DV values (should be removed)
	#print(delID <- names(temp[1,])[temp[1,]==0])
	#ID with more NA values than DV values (questionable usability)
	print(MDVoverDV <- names(temp[1,])[temp[1,]<temp[2,]])

	#Fix up MDV as it has had some errors sneak in
	datanew$MDV <- 0
	datanew$MDV[is.na(datanew$DV)] <- 1

	#Add CMT column
	datanew$CMT <- 2
	datanew$CMT[!is.na(datanew$AMT)] <- 1

	#Add OCC column
	datanew$OCC <- 1
	datanew$OCC[datanew$DAY>=4] <- 2
	datanew$OCC[datanew$DAY>=8] <- 3
	datanew$OCC[datanew$DAY>=15] <- 4

  fsex <- ifelse(datanew$GEND==1,1.23,1.04)
  datanew$CRCL <- (140-datanew$AGE)*datanew$WT*fsex/datanew$SECR  #CG w/ WT=ABW

  datanew$IBW <- ifelse(datanew$GEND==1,50+0.9*(datanew$HT-152),45.5+0.9*(datanew$HT-152))
  datanew$CRCL2 <- (140-datanew$AGE)*datanew$IBW*fsex/datanew$SECR  #CG w/ WT=IBW

  AorIBW <- with(datanew, ifelse(WT<IBW, WT, IBW))
  datanew$CRCL3 <- (140-datanew$AGE)*AorIBW*fsex/datanew$SECR  #CG w/ WT=IBW if IBW<ABW

  adjBW <- with(datanew, IBW+0.4**(WT-IBW))
  AoradjBW <- with(datanew, ifelse(WT<IBW, WT, adjBW))
  datanew$CRCL4 <- (140-datanew$AGE)*AoradjBW*fsex/datanew$SECR  #adjBW - Sawyer et. al, Leader et. al

  datanew$BMI <- datanew$WT/datanew$HT**2
  AorIoradjBW <- datanew$IBW
  AorIoradjBW[datanew$BMI < 18.5] <- datanew$WT
  AorIoradjBW[datanew$BMI < 25] <- adjBW
  datanew$CRCL5 <- (140-datanew$AGE)*AorIoradjBW*fsex/datanew$SECR  #adjBW - Winter et. al

  datanew$BLQ <- 0
  datanew$BLQ[datanew$MDV == 1 & round(datanew$TAD) == 0] <- 1
  datanew$BLQ[datanew$MDV == 1 & round(datanew$TAD) >= 20] <- 1

#Create summary tables
	datanew$DOSELVL <- dose.levels

  dataone <- lapplyBy(~X.ID, data=datanew,  oneperID)
  dataone <- bind.list(dataone)
  dim(dataone)

	dose.levelsF <- as.factor(dose.levels)
	levels(dose.levelsF) <- c("2.5mg daily","2.5mg daily (escalate to 5.0mg daily)",
														"2.5mg daily (escalate to 7.5mg daily)","5.0mg daily",
														"7.5mg daily","15mg daily","20mg daily","25mg daily",
														"30mg daily","35mg daily","50mg daily","75mg daily")
	dose.table <- data.frame("Dose Levels" = c(1:12,""),
													 "Regimen" = c(levels(dose.levelsF),"Total Patients"),
													 "No. Patients" = c(with(dataone,table(DOSELVL)),125))

	filename.out <- paste(output.dir,"dose_table.csv",sep="/")
  write.csv(dose.table, file=filename.out, row.names=F)

	temp <- with(dataone,table(GRP,DOSELVL))
	gt1 <- c(""," ",sort(unique(datanew[datanew$GRP==1,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==2,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==3,]$DOSELVL)),
					 ""," ",sort(unique(datanew[datanew$GRP==4,]$DOSELVL)),
					 ""," ",sort(unique(datanew[datanew$GRP==5,]$DOSELVL)),
					 ""," ",sort(unique(datanew[datanew$GRP==7,]$DOSELVL)),
					 " ",sort(unique(datanew[datanew$GRP==8,]$DOSELVL)))
	gt2 <- as.factor(gt1)
	levels(gt2) <- c("STUDY","Group","2.5mg daily",
									 "35mg daily","50mg daily","75mg daily",
									 "2.5mg daily (escalate to 5.0mg daily)",
									 "2.5mg daily (escalate to 7.5mg daily)",
									 "5.0mg daily","7.5mg daily","15mg daily","20mg daily",
									 "25mg daily","30mg daily")
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

### ------------------------------------- Covariate Data ------------------------------------- ###
###  			#reproducible
  file06003 <- "RAW_Clinical/datacheck_clin_06003_Output/06003_covdata.csv"
	file05115 <- "RAW_Clinical/datacheck_clin_05115_Output/05115_covdata.csv"
	file08056 <- "RAW_Clinical/datacheck_clin_08056_Output/08056_covdata.csv"
	file10016 <- "RAW_Clinical/datacheck_clin_10016_Output/10016_covdata.csv"
  data06003 <- read.csv(file06003, stringsAsFactors=F)
	data05115 <- read.csv(file05115, stringsAsFactors=F)
	data08056 <- read.csv(file08056, stringsAsFactors=F)
	data10016 <- read.csv(file10016, stringsAsFactors=F)

#FINAL datacheck and process for covspreadsheet
  datacov <- rbind(data06003,data05115,data08056,data10016)
	names(datacov)
  str(datacov)
  datacov <- orderBy(~UID+GRP, data=datacov)
	datacov$RACE[is.na(datacov$RACE)] <- 4
	datacov$RACE <- as.factor(datacov$RACE)
  levels(datacov$RACE) <- c("Caucasian","06003 Race 2","06003 Race 3","Unknown")
	datacov$DXCATNUM <- as.factor(datacov$DXCATNUM)
  levels(datacov$DXCATNUM) <- c("CLL","AML","ALL","MM")
	datacov$GEND <- as.factor(datacov$GEND)
  levels(datacov$GEND) <- c("F","M")

	filename.out <- paste(output.dir,"datacov_allstudies.csv",sep="/")
  write.csv(datacov, file=filename.out, quote=FALSE,row.names=FALSE)

#Prepare nm file
#ID TIME TAD AMT EVID OCC DV MDV ADDL II STUDY GRP DOSELVL AGE GEND WT HT SECR IBW CRCL CRCL2 CRCL3 CRCL4 CRCL5 RACE DXCAT LOQ
  nmprep <- datanew[c(1,9,10,7,8,29,12,28,13,26,27,2,4,5,15,16,17,18,24,31,30,32,33,34,35,23,21,36)]

	nmprep$WT[is.na(nmprep$WT)] <- 70
	nmprep$HT[is.na(nmprep$HT)] <- 1.75
	nmprep$HT[nmprep$HT==1.75&nmprep$GEND==0] <- 1.6
	nmprep[is.na(nmprep)] <- "."
	colnames(nmprep)[c(1,21)] <- c("#ID","RACE")

  filename.out <- "E:/Hughes/Data/PK/REDO/nmprep_allstudies.csv"
  write.csv(nmprep, file=filename.out, quote=FALSE,row.names=FALSE)

	simdata <- nmprep[!is.na(nmprep$AMT),c(1,2,3,4,5,7,9,10,11,12,13)]
	simdata[6] <- "."
	filename.out <- paste(output.dir,"simprep_allstudies.csv",sep="/")
  write.csv(simdata, file=filename.out, quote=FALSE,row.names=FALSE)
