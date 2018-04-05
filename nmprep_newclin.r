###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "nmprep_newclin"
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
  file10156 <- "PK/EXT_VAL/extval_10156.csv"
  data06003 <- read.csv(file06003, stringsAsFactors=F)
	data05115 <- read.csv(file05115, stringsAsFactors=F)
	data08056 <- read.csv(file08056, stringsAsFactors=F)
  data10016 <- read.csv(file10016, stringsAsFactors=F)
	data10156 <- read.csv(file10156, stringsAsFactors=F)

# Fix 10156 so that is binds with other data
  data10156$X.ID <- data10156$X.ID + 182
  data10156$XSAMP <- 0
  data10156$GRP <- data10156$GRP + 8
  data10156$AMT <- as.numeric(data10156$AMT)
  data10156$DV <- as.numeric(data10156$DV)
  data10156$DOSELVL <- 4
  data10156$RATE <- 0
  names(data10156)[18] <- "DXCATNUM"
  data10156$RACE <- 1
  data10156$RACE2 <- 1
  sub10156 <- data10156[-c(6, 22:23, 25, 27:29)]

#FINAL datacheck and process for nmprep
  datanew <- rbind(data06003,data05115,data08056,data10016,sub10156)
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
  dose.levels[datanew$STUDY==10156] <- 4
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

  datanew$BMI <- datanew$WT/(datanew$HT/100)**2
  AorIoradjBW <- datanew$IBW
  AorIoradjBW[datanew$BMI < 18.5] <- datanew$WT[datanew$BMI < 18.5]
  AorIoradjBW[datanew$BMI < 25] <- adjBW[datanew$BMI < 25]
  datanew$CRCL5 <- (140-datanew$AGE)*AorIoradjBW*fsex/datanew$SECR  #adjBW - Winter et. al

  datanew$BLQ <- 0
  datanew$BLQ[datanew$MDV == 1 & round(datanew$TAD) == 0] <- 1
  datanew$BLQ[datanew$MDV == 1 & round(datanew$TAD) >= 20] <- 1
  datanew$BLQ[datanew$MDV == 0 & datanew$DV < 0.00025926] <- 1
  datanew$MDV[datanew$MDV == 0 & datanew$DV < 0.00025926] <- 1

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

  dv.data <- datanew[is.na(datanew$AMT),]
  percent.na <- length(which(is.na(dv.data$DV)))/length(dv.data$DV)

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
#ID TIME TAD AMT EVID OCC DV MDV ADDL II STUDY GRP DOSELVL AGE GEND WT HT SECR IBW CRCL CRCL2 CRCL3 CRCL4 CRCL5 RACE DXCAT LOQ BSA BMI
  nmprep <- datanew[c(1,9,10,7,8,29,12,28,13,26,27,2,4,5,15,16,17,18,24,31,30,32,33,34,35,23,21,36,19,20)]
  nmprep <- rename(nmprep, c(GEND = "SEX"))

	nmprep$WT[is.na(nmprep$WT)] <- 70
	nmprep$HT[is.na(nmprep$HT)] <- 1.75
	nmprep$HT[nmprep$HT==1.75&nmprep$GEND==0] <- 1.6

  ffm.var1 <- ifelse(nmprep$SEX == 1, 6.68, 8.78)
  ffm.var2 <- ifelse(nmprep$SEX == 1, 216, 244)
  nmprep$FFM <- 9.27 * 10^3 * nmprep$WT / (ffm.var1 * 10^3 + ffm.var2 * nmprep$BMI)
	nmprep[is.na(nmprep)] <- "."
	colnames(nmprep)[c(1,26)] <- c("#ID","RACE")

  filename.out <- paste(output.dir,"nmprep_allstudies.csv",sep="/")
  write.csv(nmprep, file=filename.out, quote=FALSE,row.names=FALSE)

	simdata <- nmprep[!is.na(nmprep$AMT),c(1,2,3,4,5,7,9,10,11,12,13)]
	simdata[6] <- "."
	filename.out <- paste(output.dir,"simprep_allstudies.csv",sep="/")
  write.csv(simdata, file=filename.out, quote=FALSE,row.names=FALSE)

# Create dataset for process_auc
# Requires a dataset with times 24, 48, 72, 96 and 120 for each patient
  names(nmprep)[1] <- "ID"
  auc.capture <- ddply(nmprep, .(ID), function(x) {
    head(x, n = 5)
  })
  auc.capture$TIME <- 1:5*24
  auc.capture$TAD <- 24
  auc.capture$AMT <- "."
  auc.capture$EVID <- 0
  auc.capture$OCC <- 12
  auc.capture$DV <- "."
  auc.capture$CMT <- 2
  auc.capture$MDV <- 1
  auc.capture$ADDL <- "."
  auc.capture$II <- "."

  aucprep <- arrange(rbind(nmprep, auc.capture), ID, TIME)
  names(aucprep)[1] <- "#ID"
  filename.out <- "E:/Hughes/Data/PK/REDO/aucprep_allstudies.csv"
  write.csv(aucprep, file=filename.out, quote=FALSE,row.names=FALSE)

# Replicate the flagged dataset that was created using select_app
# These were chosen based on weighted residual and likeliness of being real
  flagprep <- nmprep
  flagprep$FLAG <- 0
  flagprep$FLAG[flagprep$ID == 6 & flagprep$OCC == 3] <- 1
  flagprep$FLAG[flagprep$ID == 6 & flagprep$OCC == 4] <- 1
  flagprep$FLAG[flagprep$ID == 7 & flagprep$OCC == 4 & flagprep$TAD < 1] <- 1
  flagprep$FLAG[flagprep$ID == 8 & flagprep$OCC == 2] <- 1
  flagprep$FLAG[flagprep$ID == 33 & flagprep$OCC == 4] <- 1
  flagprep$FLAG[flagprep$ID == 39 & flagprep$OCC == 4 & flagprep$TIME > 400] <- 1
  flagprep$FLAG[flagprep$ID == 42 & flagprep$OCC == 4 & flagprep$TAD < 3] <- 1
  flagprep$FLAG[flagprep$ID == 44 & flagprep$OCC == 4 & flagprep$TAD > 3] <- 1
  flagprep$FLAG[flagprep$ID == 48 & flagprep$OCC == 4] <- 1
  flagprep$FLAG[flagprep$ID == 48 & flagprep$OCC == 4] <- 1
  flagprep$FLAG[flagprep$ID == 50 & flagprep$OCC == 3] <- 1
  flagprep$FLAG[flagprep$ID == 63 & flagprep$OCC == 2 & flagprep$MDV == 0] <- 1
  flagprep$FLAG[flagprep$ID == 75 & flagprep$OCC == 4 & flagprep$TIME > 400] <- 1
  flagprep$FLAG[flagprep$ID == 98 & flagprep$OCC == 1 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 100 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 119 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 121 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 127 & flagprep$OCC == 1 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 129 & flagprep$OCC == 1 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 131 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 162 & flagprep$OCC == 1 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 166 & flagprep$OCC == 1 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 168 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1
  flagprep$FLAG[flagprep$ID == 173 & flagprep$OCC == 2 & flagprep$TAD == 24] <- 1

  # FIX BSA
  flagprep$BSA <- 0.007184*flagprep$WT**0.425*flagprep$HT**0.725

  filename.out <- paste(output.dir,"allnmprep_flagged_BLQ.csv",sep="/")
  write.csv(flagprep, file=filename.out, quote=FALSE,row.names=FALSE)
