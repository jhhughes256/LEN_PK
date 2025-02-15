###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_clin_08056"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)
  library(readxl)

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
  ## allinfo_excel.xls - 2 sheets of data
  #example function the sheets from this .xls were turned into .csv manually using excel (save as -> .csv)
  #allinfo.08056.data1 <- read_excel(paste(working.dir,"rawdata-lena_08056_allinfo_excel.xls",sep="/"), sheet=1)
  file.name.in <- "RAW_Clinical/rawdata-lena_08056_allinfo.csv"
  file.name.in2 <- "RAW_Clinical/rawdata-lena_08056_missingdata.csv"

  datanew <- read.csv(file.name.in, stringsAsFactors=F, na.strings=c("."))[1:12]
  datamis <- read.csv(file.name.in2,stringsAsFactors=F, na.strings=c("."))[(1:25)*2,]  #take every second value
  names(datamis) <- c("MRN","ID","DD1","CD1","DD2","CD2","DD3","CD3")

  ## demog_excel.xls - 13 sheets of data
  file.name.in3 <- "RAW_Clinical/rawdata-lena_08056_demog_excel.xlsx"
  #rows removed for double values; empty first column removed
  demog.cat <- read_excel(file.name.in3, sheet=1)[-c(5,39),-1]  #demographics
  demog.offstudy <- read_excel(file.name.in3, sheet=2)  #removed from study?
  demog.ps <- read_excel(file.name.in3, sheet=3)  #past surgery?
  demog.pt <- read_excel(file.name.in3, sheet=4)  #past treatment
  demog.bl <- read_excel(file.name.in3, sheet=5)  #toxicity around symptoms
  demog.cont <- read_excel(file.name.in3, sheet=6)  #weight, height, bsa, with accompanying cycle and visit dates
  demog.ca <- read_excel(file.name.in3, sheet=7)  #treatment information
  #demog.tox <- read_excel(file.name.in3, sheet=8)  #toxicity around bloods
  demog.haem <- read_excel(file.name.in3, sheet=9)  #haematology results - extensive
  demog.secr <- read_excel(file.name.in3, sheet=10)  #creatinine levels - extensive
  #demog.ig <- read_excel(file.name.in3, sheet=11)  #ig levels
  #demog.ul <- read_excel(file.name.in3, sheet=12)  #upper lymph nodes?
  #demog.ll <- read_excel(file.name.in3, sheet=13)  #lower lymph nodes?

  ## datapk_excel.xls - 10 sheets of data - all pk summary data
  file.name.in4 <- "RAW_Clinical/rawdata-lena_08056_datapk_excel.xlsx"
  datalena.single <- read_excel(file.name.in4, sheet=1,
    col_names=c("ID",1,"DOSE",2:10),
    col_type=c("numeric","text",rep("numeric",10)))[-c(1,2),c(1,3)]
  #datalena.dose <- read_excel(file.name.in4, sheet=2)
  #datalena.doseno <- read_excel(file.name.in4, sheet=3)
  #datalena.sum <- read_excel(file.name.in4, sheet=4)
  #datalena.sum24 <- read_excel(file.name.in4, sheet=5)
  #datalena.ave <- read_excel(file.name.in4, sheet=6)
  #dataflav.single <- read_excel(file.name.in4, sheet=7)
  #dataflav.dose <- read_excel(file.name.in4, sheet=8)
  #dataflav.sum <- read_excel(file.name.in4, sheet=9)
  #datapk.sum <- read_excel(file.name.in4, sheet=10)

  file.name.in5 <- "RAW_Clinical/rawdata-lena_08056_sigmaplot.xlsx"
  sigma.conc <- read_excel(file.name.in5, sheet=1,
    col_names = c("ID","time","conc","","ID","time","conc"),
    col_type = c(rep("numeric",7)))
  sigma.lena <- na.omit(sigma.conc[-(1:2),1:3])
  sigma.comb <- sigma.conc[-(1:2),5:7]

  #------------------------------------------------------------------------------------
#Column names ..._allinfo.csv
 #As presented
 names(datanew)

 #Sorted
 sort(names(datanew))

 #Structure
 str(datanew)

#Column names ..._missingdata.csv
 #As presented
 names(datamis)

 #Sorted
 sort(names(datamis))

 #Structure
 str(datamis)

#------------------------------------------------------------------------------------
#Observe Data

 print(temp <- with(datanew, table(PK.ID, useNA = "always")))
 #standard number of PK.ID values appears to be 40
 #some patients have 1 or 7?
 #PK.ID not currently reliable for these patients
 names(temp)[temp==7]
 names(temp)[temp==1]

 with(datanew, table(Dose.of.lena..mg., useNA = "always"))
 print(temp <- with(datanew, table(Dose.of.lena..mg.,PK.ID)))
 #Dose data only exists for these 15 patients?
 c(names(temp[1,])[temp[1,]!=0], names(temp[2,])[temp[2,]!=0])
 #Notably missing for those patients with either 1 or 7 PK.ID values
 names(temp[1,])[temp[1,]==0&temp[2,]==0]

 with(datanew, table(Treatment, useNA = "always"))
 print(temp <- with(datanew, table(Treatment,PK.ID)))
 #ID 22 has unique identifier "NOT RUN" -> explains why there is only one PK.ID value for this ID
 #These 14 patients seem standard
 names(temp[1,])[temp[1,]!=0]
 #These IDs have values for Flavo alone, but no other treatments
 names(temp[1,])[temp[1,]==0&temp[3,]!=0]
 #only have 7 values, matches with findings for PK.ID
 names(temp[1,])[temp[1,]==0&temp[3,]==7]
 #only have 1 value, matches with findings for PK.ID
 names(temp[1,])[temp[1,]==0&temp[3,]==1]


#Number of patients
 npat <- length(unique(datanew$PK.ID))-1		#negative one to account for NA
 npat
 #minus 1 from npat due to "NOT RUN" ID

#Are PK.ID on the first line of every new patient? (ASSUMPTION: PT.ID is on the first line for every new patient)
 length(unique(datanew$PT.ID))-1				#negative one to account for NA
 with(datanew, table(PK.ID,PT.ID))
 #PK.ID and PT.ID match up.. but am I missing a patient? #FIXED: account for NA
 with(datanew, table(Treatment,PT.ID, useNA = "always"))
 #PK.ID and PT.ID always appear together, even in those patients that have only 1 PK.ID value

#-------------------------------------------------------------------------------
#Convert datanew to allow for further investigation

 datanew2 <- cbind(repeat.before(datanew$PK.ID),datanew[-c(1,2)])
 colnames(datanew2) <- c("ID","TRT",colnames(datanew2[-c(1,2)]))

 #Treatment Groups
 #-----------------------------------------	n = 21
 #TRT -1	-> Flavopiridol alone				n = 20		#unique IDs c(17,18,21)
 #TRT -2	-> Flavopiridol in combination		n = 17
 #TRT 1	-> Lenalidomide alone				n = 17
 #TRT 2	-> Lenalidomide in combination		n = 17
 #TRT 0	-> "NOT RUN" (discard)				n =  1		#unique IDs c(22)
 datanew2$TRT[datanew$Treatment %in% c("Flavo alone")] <- -1
 datanew2$TRT[datanew$Treatment %in% c("Combination for Flavo")] <- -2
 datanew2$TRT[datanew$Treatment %in% c("Lena alone")] <- 1
 datanew2$TRT[datanew$Treatment %in% c("Combination for Lena")] <- 2
 datanew2$TRT[datanew$Treatment %in% c("NOT RUN")] <- 0

#Further observe the data
 #ID values
 with(datanew2, table(ID, useNA = "always"))
 with(datanew2, table(TRT, useNA = "always"))
 print(temp <- with(datanew2, table(TRT,ID)))
 #below really do only have 7 values! No more 1 values other than TRT 3
 names(temp[4,])[temp[4,]==0&temp[1,]!=0]
 #Therefore no lena data for below
 names(temp[4,])[temp[4,]==0]

#-------------------------------------------------------------------------------
#Lena subset
 datalen <- subset(datanew2,TRT>0)
 datalen$MDV <- 0
 datalen$MDV[is.na(datalen$Lena.Concentration..nM.)] <- 1

 #Check dose values
 with(datalen, table(Dose.of.lena..mg., useNA = "always"))
 print(temp <- with(datalen, table(Dose.of.lena..mg., ID)))
 #theoretically there should be 34 (17x2) dose values considering there are 17 patients each with a lena alone and lena comb regimen
 length(print(t1 <- c(names(temp[1,])[temp[1,]==2],names(temp[2,])[temp[2,]==2])))
 length(print(t2 <- c(names(temp[1,])[temp[1,]==1],names(temp[2,])[temp[2,]==1])))
 #only 4 patients have 2 values and 8 patients have 1 value -> 12/17 patients with 16/34 data points
 #IDs with dose values
 len.dose <- sort(as.numeric(c(t1,t2)))

 print(temp <- with(datalen, table(is.na(Lena.Concentration..nM.),ID,TRT)))
 #Lena alone - see below			- 14/17 patients
 length(print(len.alon <- c(names(temp[1,,1])[temp[1,,1]!=0])))
 #Lena comb  - see below			- 14/17 patients
 length(print(len.comb <- c(names(temp[1,,2])[temp[1,,2]!=0])))

 print(temp <- with(datalen, table(MDV,ID)))
 #14/17 patients with concentration data (more than dose data)				#identical to flav data
 length(print(len.conc <- as.numeric(names(temp[1,])[temp[1,]!=0])))
 #patients with concentration values but no dose values
 len.conc[!len.conc %in% len.dose]
 #patients with dose values but no concentration values 					#identical vector to flav data
 len.dose[!len.dose %in% len.conc]
 #Lena alone - missing some values on patients 9 (1 MDV) and 19 (7 MDV)
 with(datalen, table(MDV,ID,TRT))

 print(temp <- with(datalen, table(Lena.Time..hr., ID, useNA = "always")))
 #no time data for patients below -> this is ok there is no dose or conc data for these patients
 print(lena.time <- as.numeric(names(temp[dim(temp)[1],])[temp[dim(temp)[1],]!=0]))

 print(temp <- with(datalen, table(TRT, Lena.Time..hr., useNA = "always")))
 #TRT 1 time points
 names(temp[1,])[temp[1,]!=0][!is.na(names(temp[1,])[temp[1,]!=0])]
 #TRT 2 time points
 names(temp[2,])[temp[2,]!=0][!is.na(names(temp[2,])[temp[2,]!=0])]

 #Found when checking "empty" columns
 with(datalen, table(Flavo.Time..hr., useNA = "always"))
 #problematic, single time + conc datapoint from Flavo is included in subset ID15

#Need to separate AMT values from OBS values for sake of MDV, EVID etc.
 datalen$TRT <- as.numeric(datalen$TRT)
 temp <- orderBy(~TRT+ID+Flavo.Time..hr., data=datalen)
 temp <- lapplyBy(~ID, data=temp, oneperID)
 t1 <- bind.list(temp)
 temp <- orderBy(~-TRT+ID+Flavo.Time..hr., data=datalen)
 temp <- lapplyBy(~ID, data=temp, oneperID)
 t2 <- bind.list(temp)

 #remove empty dose row
 len.amt <- rbind(t1[!is.na(t2[11]),])
 len.amt$Lena.Concentration..nM. <- NA
 len.amt$MDV <- 1

 #add in correct dosage from PK summary
 lendoses <- lapplyBy(~ID, data=datalena.single[datalena.single$ID<=25,],oneperID)
 lendoses <- na.omit(bind.list(lendoses))
 len.amt$Dose.of.lena..mg. <- lendoses$DOSE


#Setup datalen for rbind with len.amt
#Needs Rate values <- 0 (as dataframe only contains observations)
#Needs Dose values <- NA (	""		""		""		""	)
 datalen$Dose.of.lena..mg. <- NA
 datalen$Lena.Dose..mmol. <- NA
 datalen <- rbind(datalen,len.amt)
 datalen$RATE <- 0

#-------------------------------------------------------------------------------
#Flav subset
 datafla <- subset(datanew2,TRT<0)
 datafla$MDV <- 0
 datafla$MDV[is.na(datafla$Flavo.Concentration..nM.)] <- 1

 #Check dose values
 with(datafla, table(total.Flavo..mg., useNA = "always"))
 print(temp <- with(subset(datafla,!is.na(Flavo.Concentration..nM.)), table(TRT,ID)))
 #Flavo alone - see below			- 20/20 patients
 length(print(t1 <- names(temp[1,][temp[1,]!=0])))
 #Flavo comb  - see below			- 14/17 patients
 length(print(t2 <- names(temp[2,][temp[2,]!=0])))
 #IDs with dose values
 flav.dose <- sort(unique(as.numeric(c(t1,t2))))
 #there should be 37 dose values
 #only 3 patients have 2 values and 9 patients have 1 value -> 12/20 patients with 15/37 data points
 length(t1[t1 %in% t2])

 print(temp <- with(datafla, table(MDV,ID)))
 #all 20 patients have concentration data for flavopiridol
 length(print(as.numeric(names(temp[1,])[temp[1,]!=0])))
 print(temp <- with(datafla, table(MDV,ID,TRT)))

 #flavo comb patients with concentration data		- 14/17 patients
 length(print(flav.conc <- as.numeric(names(temp[1,,2])[temp[1,,2]!=0])))
 #identical to vector for len data
 flav.conc[!flav.conc %in% flav.conc]
 #patients with concentration values but no dose values
 flav.conc[!flav.conc %in% as.numeric(t2)]
 #patients with dose values but no concentration values 					#identical vector to flav data
 as.numeric(t2)[!as.numeric(t2) %in% flav.conc]

 print(temp <- with(subset(datafla,TRT==-1), table(Flavo.Time..hr., ID, useNA = "always")))
 #time values for all flav alone patients except 1x48hr time -> explained by lena subset
 names(temp[dim(temp)[1],])[temp[dim(temp)[1],]!=0]

 print(temp <- with(subset(datafla,TRT==-2), table(Flavo.Time..hr., ID, useNA = "always")))
 #no time values for ID below
 print(flav.time <- as.numeric(names(temp[dim(temp)[1],])[temp[dim(temp)[1],]!=0]))
 #identical to lena data
 flav.time[!flav.time %in% lena.time]

 #no Lena data in this subset

#Alter flav dataset to fix missing 48hr value that was located in lena dataset
 temp <- subset(datalen,ID==15)
 #row 420 contains the stray values
 temp[6]==48
 temp <- subset(datafla,ID==15)
 #must replace row 401 in datanew
 is.na(temp[6])
 #this equates to row 191 in datafla
 which(rownames(datafla)==401,)

 #WARNING: Not ordered by TIME
 #Added zero to end of datanew2 as it does not include MDV
 #col 2 of datafla changed as its TRT is set to that of lena
 datafla[191,] <- c(datanew2[420,],0)
 datafla[191,2] <- -1

#Alter flav dataset to have AMT at correct time (first dose at TIME=0, second dose at TIME=0.5
#t1 is dataset containing AMT values for each ID for Flav alone
#t2 is		""		""		""		for Flav comb
 datafla$TRT <- as.numeric(datafla$TRT)
 temp <- orderBy(~-TRT+ID+Flavo.Time..hr., data=datafla)
 temp <- lapplyBy(~ID, data=temp, oneperID)
 t1 <- bind.list(temp)
 temp <- orderBy(~TRT+ID+Flavo.Time..hr., data=datafla)
 temp <- lapplyBy(~ID, data=temp, oneperID)
 t2 <- bind.list(temp)[1:17,]			#final 3 rows excluded as only 17 patients for flav comb
 temp <- orderBy(~-TRT+ID+Flavo.Time..hr., data=datalen)
 temp <- lapplyBy(~ID, data=temp, oneperID)
 t3 <- bind.list(temp)

 #remove empty dose row
 t1 <- t1
 t2 <- t2[!is.na(t3[11]),]

 t1$Dose.of.Flavo..mg.[t1$Dose.of.Flavo..mg.==""] <- "NA/NA"
 t2$Dose.of.Flavo..mg.[t2$Dose.of.Flavo..mg.==""] <- "NA/NA"

 #split dose values into pairs then select second value of the pair
 dose1 <- as.numeric(unlist(strsplit(t1$Dose.of.Flavo..mg.,"/")))
 dose2 <- as.numeric(unlist(strsplit(t2$Dose.of.Flavo..mg.,"/")))

 #Create AMT dataframe to rbind to datafla
 flav.amt <- data.frame("ID" 						= c(rep(t1$ID,each=2),rep(t2$ID,each=2)),
            "TRT"						= c(rep(t1$TRT,each=2),rep(t2$TRT,each=2)),
            "Dose.of.Flavo..mg." 		= c(dose1,dose2),
            "total.Flavo..mg." 		= c(as.vector(rbind(t1$total.Flavo..mg.,rep(0,length(t1$ID)))),as.vector(rbind(t2$total.Flavo..mg.,rep(0,length(t2$ID))))),
            "Total.Flavo..mmol."		= c(as.vector(rbind(t1$Total.Flavo..mmol.,rep(0,length(t1$ID)))),as.vector(rbind(t2$Total.Flavo..mmol.,rep(0,length(t2$ID))))),
            "Flavo.Time..hr."			= c(rep(c(0,0.5),times=length(t1$ID)),rep(c(0,0.5),times=length(t2$ID))),
            "Flavo.Concentration..nM."	= NA,
            "Dose.of.lena..mg."		= NA,
            "Lena.Dose..mmol."			= NA,
            "Lena.Time..hr."			= NA,
            "Lena.Concentration..nM."	= NA,
            "MDV"						= 1,
            "RATE"						= c(as.vector(rbind(rep(0,length(t1$ID)),dose1[1:(length(dose1)/2)*2]/4)),as.vector(rbind(rep(0,length(t2$ID)),dose1[1:(length(dose2)/2)*2]/4))))

#Setup datafla for rbind with flav.amt
#Needs Rate values <- 0 (as dataframe only contains observations)
#Needs Dose values <- NA (	""		""		""		""	)
 datafla$RATE <- 0
 datafla$Dose.of.Flavo..mg. <- NA
 datafla$total.Flavo..mg. <- NA
 datafla$Total.Flavo..mmol. <- NA

 datafla <- rbind(datafla,flav.amt)

#------------------------------------------------------------------------------
#Setup covariate data; MRN's required for matching
  IDcor <- datamis[1:2]
  demogcat2 <- merge(IDcor,demog.cat)[c(1:2,9:11)]
  names(demogcat2)[3] <- "Age"

  names(demog.cont)[c(1,3,7,9,11)] <- c("MRN","CYCLE","WT","HT","BSA")
  demogcont1 <- merge(IDcor,demog.cont)[-c(3,5:7,9,11,13)]
  print(temp <- with(demogcont1, table(ID, CYCLE, useNA="always")))
  print(temp2 <- names(temp[,2][temp[,2]==0]))  #no values for cycle 2 for some individuals
  len.conc[which(len.conc %in% temp2)]    #does not match with len.conc, not individuals that have data
  names(temp[,9][temp[,9]!=0])  #individuals with NAs
  demogcont2 <- demogcont1[demogcont1$ID %in% as.numeric(names(temp[,2][temp[,2]!=0])),]    #remove from data
  demogcont2 <- demogcont2[!is.na(demogcont2$WT)&!is.na(demogcont2$CYCLE),]

  names(demog.secr)[1] <- "MRN"
  demogsecr1 <- merge(IDcor,demog.secr)
  with(demogsecr1, table(Day, Cycle, useNA="always"))
  #all of cycle 1 is considered to be day 1 despite having multiple different visit dates
  #have multiple day 1 and day 3 for cycle 2 (lena cycle)
  #I only care about cycle 2 day 1 and day 3
  with(na.omit(demogsecr1[demogsecr1$Cycle==2,]), table(Day, ID, useNA="always"))
  #patients have secr on at least either day 1 or day 3

  demogsecr1$Value <- as.numeric(demogsecr1$Value)
  demogsecr1$Day <- as.numeric(demogsecr1$Day)
  demogsecr2 <- na.omit(demogsecr1[demogsecr1$Cycle==2&demogsecr1$Day<=3,])

  demogsecr3 <- ddply(demogsecr2, .(ID), function(x){
    df1 <- x[x$Day==1,]
    v1 <- ifelse(length(df1$Value)!=0,mean(df1$Value),NA)
    df3 <- x[x$Day==3,]
    v3 <- ifelse(length(df3$Value)!=0,mean(df3$Value),NA)
    if(length(df1$Value)==0){
      df <- rbind(df3[rep(1,3),1:6])
    }else{
      if(length(df3$Value)==0){
        df <- rbind(df1[rep(1,3),1:6])
      }else{
        df <- rbind(df1[c(1,1),1:6],df3[1,1:6])
      }
    }
    temp <- c(v1,v1,v3)
    secr <- ifelse(length(temp[!is.na(temp)])==2,temp[1],temp[!is.na(temp)])
    data.frame(df,"DAY" = c(1,2,3),"SECR" = secr)
  })[c(1:2,7:8)]

  datalen1 <- merge(IDcor,datalen)
  datalen2 <- merge(datalen1,demogcat2)

  datafla1 <- merge(IDcor,datafla)
  datafla2 <- merge(datafla1,demogcat2)

  sigma.id <- unique(sigma.lena$ID)[!unique(sigma.lena$ID) %in% unique(as.numeric(datalen2$ID))]
  #sigma.lena1 <-
  #sigma.comb1 <-

#------------------------------------------------------------------------------
  #Convert datanew2 to old format
  #Additionally stack flav and lena AMT, TIME, DV

  #Treatment Groups
  #-----------------------------------------	n = 21
  #TRT -1	-> Flavopiridol alone				n = 20		#unique IDs c(17,18,21)
  #TRT -2	-> Flavopiridol in combination		n = 17
  #TRT 1	-> Lenalidomide alone				n = 17
  #TRT 2	-> Lenalidomide in combination		n = 17
  #TRT 0	-> "NOT RUN" (discard)				n =  1		#unique IDs c(22)

  dataint <- data.frame(
    "ID" = c(datafla2$ID,datalen2$ID),
    "MRN" = c(datafla2$MRN,datalen2$MRN),
    "STUDY" = 8056,
    "TRT" = c(datafla2$TRT,datalen2$TRT))

  dataint$AMT <- c(datafla2$Dose.of.Flavo..mg.,datalen2$Dose.of.lena..mg.)

  dataint$RATE <- c(datafla2$RATE,datalen2$RATE)

  dataint$TIME <- c(datafla2$Flavo.Time..hr.,datalen2$Lena.Time..hr.)

  dataint$DV <- c(datafla2$Flavo.Concentration..nM.*401.84,datalen2$Lena.Concentration..nM.*259.26)		#molar mass for flavo and lena
  dataint$DV[dataint$DV<=0] <- NA
	dataint$DV <- dataint$DV/(1000*1000)																																		#conversion from ng/L to ug/mL

  dataint$MDV <- c(datafla2$MDV,datalen2$MDV)
  dataint$MDV[dataint$DV<0.00025926] <- 1

  dataint$LNDV <- log(dataint$DV)

  dataint$DAY <- 3
  dataint$DAY[dataint$TRT==-1] <- 1
  dataint$DAY[dataint$TRT==1]  <- 2

  dataint$CYCLE <- 2
  dataint$CYCLE[dataint$TRT==-1] <- 1

  dataint$DOSEMG <- c(datafla2$total.Flavo..mg.,datalen2$Dose.of.lena..mg.)
  dataint$DOSEMG[dataint$DOSEMG==0] <- NA

  #dataint$BLQ <- 0
  #dataint$BLQ[datanew$NOTE == "DV_BLQ"] <- 1
  #with(dataint, table(BLQ, useNA = "always"))
  #dataint$DV[dataint$BLQ==1] <- NA

  #dataint$DNUM <- datanew$Dose.no
  #dataint$DNUM <- unlist(lapplyBy(~ID, data=dataint, function(d) impute(d$DNUM)))

  #dataint$OCC <- dataint$DNUM

  dataint$AGE <- c(datafla2$Age,datalen2$Age)

  dataint$GEND <- c(datafla2$Gender,datalen2$Gender)				  	#1 is male, 0 is female
  dataint$GEND[dataint$GEND %in% c("M")] <- 1
  dataint$GEND[dataint$GEND %in% c("F")] <- 0

  #DXCATNUM == 1 -> (Relapsed/Refactory) Chronic Lymphocytic Leukaemia - Maddocks K. et al. 2015
  dataint$DXCATNUM <- 1

  #White              	1
	#African American			2

  dataint$RACE2 <- c(datafla2$Race,datalen2$Race)
  dataint$RACE2[dataint$RACE2 %in% c("White")] <- 1
  dataint$RACE2[dataint$RACE2 %in% c("Black or African American")] <- 2

  datanew3 <- merge(dataint,demogcont2)
  datanew3 <- merge(datanew3,demogsecr3)
  datanew3$ID <- as.numeric(levels(datanew3$ID))[datanew3$ID]

 #-----------------------------------------------------------------



 #-----------------------------------------------------------------
  dataall <- datanew3

  dataall <- orderBy(~ID+TRT+TIME+AMT, data=dataall)

  len.alon <- as.numeric(len.alon)
  len.comb <- as.numeric(len.comb)

  dataall.lena <- rbind(subset(dataall,TRT==1)[subset(dataall,TRT==1)$ID %in% len.alon,],subset(dataall,TRT==2)[subset(dataall,TRT==2)$ID %in% len.comb,])
  dataall.lena$DOSEMG <- repeat.before(dataall.lena$DOSEMG)

#-------------------------------------------------------------------------------
# Check subject numbers
   with(dataall.lena, table(ID))

#----------------------------------------------------------------------------------------------------------------------
#Calculate dose normalized concentrations

   #Check distribution of DV
    plotobj <- qplot(x=DV, geom="histogram", data=dataall.lena)
    plotobj

    filename.out <- paste(output.dir,"Histogram_DV",sep="/")
    suppressWarnings(to.png(plotobj,filename.out))
    #DV has a wide range of values - 4 orders of magnitude at least

   #Check distribution of log DV
    plotobj <- qplot(x=LNDV, geom="histogram", data=dataall.lena)
    plotobj

    filename.out <- paste(output.dir,"Histogram_DVlog",sep="/")
    suppressWarnings(to.png(plotobj,filename.out))
    #DV has a wide range of values - 4 orders of magnitude at least


   #Calculate dose normalised DV
   #Units are ng/ml per mg

    dataall.lena$DVNORM <- dataall.lena$DV/dataall.lena$DOSEMG

#----------------------------------------------------------------------------------------------------------------------
#Count missing covariate data
#Missing by Study

 covnames <- as.formula("~AGE+GEND+WT+HT+BSA+RACE2+DXCATNUM+SECR")
 covdata <- subset(dataall.lena, select=c("ID","DOSEMG","AGE","GEND","WT","HT","BSA","RACE2","DXCATNUM","SECR"))

#Reassign missing
 covdata[covdata==-1] <- NA

  #finish off
 missingbystudy <- ddply(covdata, .(DOSEMG), colwise(calculate.percent.missing))

 filename.out <- paste(output.dir,"Missing_by_Group.csv",sep="/")
 write.csv(missingbystudy, file=filename.out, row.names=F)



#Missing by Subject
 missingbysubject <- ddply(covdata, .(DOSEMG,ID), colwise(calculate.percent.missing))
 filename.out <- paste(output.dir,"Missing_by_Subject.csv",sep="/")
 write.csv(missingbysubject, file=filename.out, row.names=F)


#-------------------------------------------------------------------------------
#Subset covariates

 ##Keeps missing as -1 - use for categorical summary
 dataallone <- lapplyBy(~ID, data=dataall.lena,oneperID)
 dataallone <- bind.list(dataallone)
 #dim(dataallone)

 #with(dataallone,table(RATE,useNA="always"))View

 ##Sets missing to NA - use for continuous summary
 covdataone <- lapplyBy(~ID, data=covdata,  oneperID)
 covdataone <- bind.list(covdataone)
 #dim(covdataone)


 dataallone$RACEf <- factor(dataallone$RACE2, labels=c("White or Caucasian","Other"))

 dataallone$SEXf <- as.factor(dataallone$GEND)
 levels(dataallone$SEXf) <- c("female","male")

 dataallone$IDf <- as.factor(dataallone$ID)

 dataallone$DXCAT2f <- as.factor(dataallone$DXCATNUM)
 levels(dataallone$DXCAT2f) <- "CLL"

 #DXCAT2 "Relapsed Chronic Lymphocytic Leukaemia" <- 1
 #DXCAT2 "Acute Myeloid Leukaemia" <- 2
 #DXCAT2 "Acute Lymphoblastic Leukaemia" <- 3

#-----------------------------------------------------------------------
#Summary of study characteristics

 #Do all subjects have PK data
 DVtest <- summaryBy(DV ~ ID, data=dataall.lena, FUN=mean, na.rm=T)
 DVtest <- DVtest[is.na(DVtest$DV.mean)==T,]
 DVtestID <- DVtest$ID
 DVtestID
 #All patients have PK data


 #Do all subjects have dose data
 AMTtest <- summaryBy(AMT ~ ID, data=dataall.lena, FUN=mean, na.rm=T)
 AMTtest <- AMTtest[is.na(AMTtest$AMT.mean)==T,]
 AMTtestID <- AMTtest$ID
 AMTtestID
 with(dataallone,table(AMT,useNA="always"))
  #All subjects have dose data

 #Do all subjects have nmRATE data
 #RATEtest <- summaryBy(RATE ~ ID, data=dataall, FUN=mean, na.rm=T)
 #RATEtest <- RATEtest[is.na(RATEtest$RATE.mean)==T,]
 #RATEtestID <- RATEtest$ID
 #RATEtestID
 #with(dataallone,table(RATE,useNA="always"))
  #All subjects have nmRATE data but not in the right place in some cases

 #DV count by Study
 #Calculates data for Report Table 1
 DVcount <- summaryBy(DV ~ TRT, data=dataall.lena, FUN=lengthNA)
 names(DVcount) <- c("Group","DVcount")
 DVcount


 #Subject count by Study
 SUBcount <- ddply(dataall.lena, .(TRT), function(df) count.unique(df$ID))
 names(SUBcount) <- c("Group","SUBcount")
 SUBcount

 #Dose count by Study
 AMTcount <- ddply(dataall.lena, .(TRT), function(df) lengthNA(df$AMT))
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
 DVcount <- summaryBy(DV ~ TRT+DOSEMG, data=dataall.lena, FUN=lengthNA)
 names(DVcount) <- c("Group","Dose Level","DVcount")
 DVcount


#Subject count by Group and DoseLevel
 SUBcount <- ddply(dataall.lena, .(TRT,DOSEMG), function(df) count.unique(df$ID))
 names(SUBcount) <- c("Group","Dose Level","SUBcount")
 SUBcount


#Dose count by Group and DoseLevel
 AMTcount <- ddply(dataall.lena, .(TRT,DOSEMG), function(df) lengthNA(df$AMT))
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
  DVcountdata <-  ddply(dataall.lena, .(TRT,ID,DAY), function(df) DV.present(df$DV))


  withDVbyGRPWEEK <- ddply(DVcountdata, .(TRT,DAY), function(df) sum(df$V1))  #GOLD
  withDVbyGRPWEEK


  #Not all subjects with Day 1 data also have Day 4 data
  filename.out <- paste(output.dir,"DVwith_group_week.csv",sep="/")
  write.csv(withDVbyGRPWEEK, file=filename.out)

#----------------------------------------------------------------------------------------------------------------------
#Subset some plot data

  plotdata <- subset(dataall.lena)
  BINnumber <- 3

  plotdata$DOSEMGf <- as.factor(plotdata$DOSEMG)

  plotdata$TRTf <- as.factor(plotdata$TRT)
  levels(plotdata$TRTf) <- paste("Group",levels(plotdata$TRTf))

 #plotdata$VOSf <- as.factor(plotdata$VOS)
 #levels(plotdata$VOSf) <- c("Placebo","Adrug")

  plotdata$DAYf <- as.factor(plotdata$DAY)
  levels(plotdata$DAYf) <- paste("Week",levels(plotdata$DAYf))

 plotdata$SEXf <- as.factor(plotdata$GEND)
 levels(plotdata$SEXf) <- c("female","male")

 plotdata$RACEf <- as.factor(plotdata$RACE2)
 levels(plotdata$RACEf) <- c("White or Caucasian","Other")

 plotdata$DOSELVLf <- as.factor(plotdata$DOSEMG)
 levels(plotdata$DOSELVLf) <- paste("Dose Level",levels(plotdata$DOSELVLf))

 #plotdata$DOSE_bin <- cut2(plotdata$DOSEMG, g=BINnumber)

 plotdata$AGE_bin <- cut2(as.numeric(plotdata$AGE), g=BINnumber)

 plotdata$WT_bin <- cut2(as.numeric(plotdata$WT), g=BINnumber)

 plotdata$HT_bin <- cut2(as.numeric(plotdata$HT), g=BINnumber)

 plotdata$BSA_bin <- cut2(as.numeric(plotdata$BSA), g=BINnumber)

 #plotdata$BMI_bin <- cut2(plotdata$BMI, g=BINnumber)

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
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSEMGf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose")
  plotobj

  filename.out <- paste(output.dir,"Overview_ConcObs_vs_TIME_by_DOSE",sep="/")
  to.png(plotobj,filename.out)

 #Conc vs TIME Week 1
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=TRTf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Treatment")
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_TRT",sep="/")
  to.png(plotobj,filename.out)

 #Conc vs TIME Week 1 per ID
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=TRTf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Treatment")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"Week1_ConcObs_vs_TIME_by_ID",sep="/")
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
  plotByFactor("TRTf","Treatment Group")
  plotByFactor("SEXf","SEX")
  plotByFactor("DOSEMGf","Dose Level")
  #plotByFactor("WEEKf","Week")
  #plotByFactor("DOSE_bin","Binned Dose (mg)")
  plotByFactor("AGE_bin","Binned Age (years)")
  plotByFactor("WT_bin","Binned Weight (kg)")
  plotByFactor("HT_bin","Binned Height (kg)")
  plotByFactor("BSA_bin","Binned BSA (m2)")
  #plotByFactor("BMI_bin","Binned BMI (kg per m2)")
  plotByFactor("DXCAT2f","Disease category")


#-------------------------------------------------------------------------
#Summarize Categorical covariates

 #dataallone$SUMCOL <- "All Groups"

 #Sex
 #covCatSexStudy <- with(dataallone,ftable(GRP,SEXf, useNA="ifany", dnn=c("GRP","CATEGORY")))
 #covCatSexStudy  <- data.frame("COV"="SEX",covCatSexStudy)


 ##Race
 #covCatRaceStudy <- with(dataallone,ftable(GRP,RACEf, useNA="ifany", dnn=c("GRP","CATEGORY")))
 #covCatRaceStudy  <- data.frame("COV"="RACE",covCatRaceStudy)


 ##DXcategory
 #covCatDXcatStudy <- with(dataallone,ftable(GRP,DXCAT2f, useNA="ifany", dnn=c("GRP","CATEGORY")))
 #covCatDXcatStudy <- data.frame("COV"="DXCAT",covCatDXcatStudy)


 ##Collate
 #covCatTable <- rbind(covCatSexStudy,covCatRaceStudy,covCatDXcatStudy)
 #covCatTable

 ##Return to original order
 #covCatTable <- orderBy(~COV, covCatTable)  #GOLD - sort by factor levels

    #Define reassignment of column names (if any) rtf allowed
#    colNames <- c(
#                  "COV","Covariate Code",
#                  "CATEGORY","Category",
#                   "RACE","Race",
#                   "SEX","Gender",
#                   "Freq","Count",
#				  "DXCAT","Disease"
#                  )
#
 #colData <- data.frame(matrix(colNames, byrow=T, ncol=2),stringsAsFactors=F)
 #names(colData) <- c("ColNameOld","ColNameNew")

   #Reassign column names
 #covCatTable$COV <- gsub.all(colData$ColNameOld,colData$ColNameNew,covCatTable$COV)
 #names(covCatTable) <- gsub.all(colData$ColNameOld,colData$ColNameNew,names(covCatTable))

   #Write as an rtf file
 #filename.out <- paste(output.dir,"CatCovSummary",sep="/")
 #write.csv(covCatTable, file=filename.out)

#-----------------------------------------------------------------------------------------------------
   #Pairs-plot of continuous covariates
 #plotobj <- NULL
 #covdataonecont <- subset(covdataone, select=c("AGE","WT","BSA","BMI","SECR"))
 #plotobj <- ggpairs(na.omit(covdataonecont))

 #filename.out <- paste(output.dir,"Overview_cont_cov_pairs",sep="/")
 #to.png(plotobj,filename.out)


#-----------------------------------------------------------------------------------------------------
 #Pairs-plot of categorical covariates

 #covcatdataonecat <- subset(dataallone, select=c("SEXf","RACEf","DXCAT2f"))

 #plotobj <- NULL
 #plotobj <- ggpairs(na.omit(covcatdataonecat))

 #filename.out <- paste(output.dir,"Overview_cat_cat_pairs",sep="/")
 #to.png(plotobj,filename.out)


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
# theme_bw2 <- theme_set(theme_bw(base_size = 20))
# theme_bw2 <- theme_update(plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
# axis.title.x=element_text(size = 18, vjust = 0),
# axis.title.y=element_text(size = 18, vjust = 1, angle = 90),
# strip.text.x=element_text(size = 14),
# strip.text.y=element_text(size = 14, angle = 90))


plotIndexCont <- function(CovColname,CovText)
{
  #Debug
  #CovColname <- "HT"
  #CovText <- "Height (cm)"

  plotobj <- ggplot(data=dataallone)
  plotobj <- plotobj + geom_point(aes_string(y=as.numeric(CovColname), x="IDf"), size=3)
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
#plotIndexCont("AGE","Age~(years)")
#plotIndexCont("WT","Weight~(kg)")
#plotIndexCont("BSA","Body~Surface~Area~(m^2)")
#plotIndexCont("BMI","Body~Mass~Index~(kg/m^2)")
plotIndexCat("GEND","Patient~Sex")
plotIndexCat("RACE2","Patient~Race")
plotIndexCat("DXCATNUM","Diagnosis~Category")
plotIndexCont("SECR","Serum~Creatinine~(umol/L)")

#Data prep
	# [1] "#ID"      "STUDY"    "XSAMP"    "GRP"      "DOSELVL"  "DOSEMG"   "AMT"      "RATE"     "TIME"
	#[10] "TAD"      "DAY"      "DV"       "MDV"      "LNDV"     "AGE"      "GEND"     "WT"       "HT"
	#[19] "BSA"      "BMI"      "DXCATNUM" "RACE"     "RACE2"    "SECR"     "DVNORM"   "ADDL"     "II"

	dataFIX <- data.frame("ID" = (dataall.lena$ID+96), "STUDY" = dataall.lena$STUDY)

  dataFIX$XSAMP <- 0
  dataFIX$GRP <- 5
	dataFIX$GRP[dataall.lena$TRT==2] <- 6

  dataFIX$DOSELVL <- 1
	dataFIX$DOSELVL[dataall.lena$DOSEMG==5] <- 2

  dataFIX$DOSEMG <- as.numeric(dataall.lena$DOSEMG)
	dataFIX$AMT <- dataall.lena$AMT
  dataFIX$RATE <- dataall.lena$RATE
  dataFIX$TIME <- dataall.lena$TIME+24
	dataFIX$TIME[!is.na(dataFIX$AMT)] <- 0
	dataFIX$TAD <- dataall.lena$TIME
	dataFIX$TAD[dataFIX$TAD==48] <- 24

  dataFIX$DAY <- dataall.lena$DAY
	dataFIX$DAY[dataFIX$TIME==0] <- 1
	dataFIX$TIME[dataFIX$DAY==3] <- dataFIX$TIME[dataFIX$DAY==3]+24
	dataFIX$DAY[dataFIX$TIME==48] <- 2
	dataFIX$DAY[dataFIX$TIME==96] <- 4

  dataFIX$DV <- dataall.lena$DV
  dataFIX$MDV <- dataall.lena$MDV
  dataFIX$LNDV <- dataall.lena$LNDV

  dataFIX$AGE <- dataall.lena$AGE
  dataFIX$GEND <- dataall.lena$GEND
  dataFIX$WT <- dataall.lena$WT
  dataFIX$HT <- dataall.lena$HT
  dataFIX$BSA <- dataall.lena$BSA
  dataFIX$BMI <- as.numeric(dataall.lena$WT)/as.numeric(dataall.lena$HT)**2
  dataFIX$DXCATNUM <- dataall.lena$DXCATNUM
  dataFIX$RACE <- dataall.lena$RACE2
  dataFIX$RACE2 <- dataall.lena$RACE2
  dataFIX$SECR <- dataall.lena$SECR*88.4	#convert from mg/dL to umol/L
  dataFIX$DVNORM <- dataall.lena$DVNORM

  dataFIX$ADDL <- NA
  dataFIX$ADDL[!is.na(dataFIX$AMT)] <- 20

  dataFIX$II <- NA
	dataFIX$II[!is.na(dataFIX$AMT)] <- 24

	dataFIX <- orderBy(~ID+TIME+GRP+AMT, data=dataFIX)
	colnames(dataFIX)[1] <- "#ID"

	filename.out <- paste(output.dir,"08056_finaldata.csv",sep="/")
  write.csv(dataFIX, file=filename.out, row.names=FALSE)

#------------------
#Covariate data
	# [1] "UID"      "ID"       "STUDY"    "GRP"      "DOSELVL"  "DOSEMG"   "AGE"      "GEND"
	# [8] "WEIGHTLB" "HEIGHTFT" "DXCATNUM" "RACE"     "SECRMGDL"

	dataCOV <- data.frame("UID" = dataallone$ID+96, "ID" = dataallone$ID, "STUDY" = dataallone$STUDY)

  dataCOV$GRP <- dataallone$TRT+4
  dataCOV$DOSELVL <- 1
  dataCOV$DOSEMG <- dataallone$AMT
	dataCOV$DOSELVL[dataCOV$DOSEMG==5] <- 2
  dataCOV$AGE <- dataallone$AGE
  dataCOV$GEND <- dataallone$GEND
  dataCOV$WT <- as.numeric(dataallone$WT)
  dataCOV$HT <- as.numeric(dataallone$HT)
  dataCOV$DXCATNUM <- dataallone$DXCATNUM
  dataCOV$RACE <- dataallone$RACE2
  dataCOV$SECRMGDL <- dataallone$SECR
	dataCOV$MDV <- dataallone$MDV
	dataCOV$DOSEMG[dataCOV$MDV==0&is.na(dataCOV$DOSEMG)] <- "Missing"

	dataCOV <- dataCOV[!is.na(dataCOV$DOSEMG),]

	filename.out <- paste(output.dir,"08056_covdata.csv",sep="/")
  write.csv(dataCOV[-dim(dataCOV)[2]], file=filename.out, row.names=FALSE)


# ------------------------------------------------------------------------------
# Clean Data (non-nmprep)
  library(dplyr)
  clean.data <- dataFIX
  names(clean.data)[1] <- "ID"
  clean.data <- clean.data %>%
    select(c(ID, STUDY, DOSEMG, TIME, TAD, DAY, DV, MDV,
      AGE, GEND, WT, HT, DXCATNUM, RACE2, SECR)) %>%
    rename(DXCAT = DXCATNUM, SEX = GEND, DVMGL = DV, WTKG = WT, HTCM = HT, SECRUMOLL = SECR, RACE = RACE2)

  clean.data$ID <- clean.data$ID - 96

  dxcat.l <- c("CLL", "AML", "ALL")
  clean.data$DXCAT <- factor(clean.data$DXCAT, levels = c(1, 2, 3))
  levels(clean.data$DXCAT) <- dxcat.l

  sex.l <- c("F", "M")
  clean.data$SEX <- factor(clean.data$SEX, levels = c(0, 1))
  levels(clean.data$SEX) <- sex.l

  race.l <- c("W", "B")
  clean.data$RACE <- factor(clean.data$RACE, levels = c(1, 2))
  levels(clean.data$RACE) <- race.l

  filename.out <- paste(output.dir,"08056_cleandata.csv",sep="/")
  write.csv(clean.data, file=filename.out, row.names=FALSE)
