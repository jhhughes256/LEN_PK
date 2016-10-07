###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_clin_10156"
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
  file.name.in1 <- "RAW_Clinical/rawdata-lena_10156_datanew_excel.xlsx"
  datanew.pkd <- read_excel(file.name.in1, sheet=1)  #pk data
  #datanew.avg <- read_excel(file.name.in1, sheet=2)  #average cp
  #datanew.plot <- read_excel(file.name.in1, sheet=3)  #plots
  #datanew.relt <- read_excel(file.name.in1, sheet=4)  #pk data relative time
  #datanew.avre <- read_excel(file.name.in1, sheet=5)  #pk data relative time (looks like this contains errors)
  datanew.act <- read_excel(file.name.in1, sheet=6)  #actual dose time

  file.name.in2 <- "RAW_Clinical/rawdata-lena_10156_datapk_excel.xlsx"
  datapk.demog <- read_excel(file.name.in2, sheet=1)  #demographic data
  #datapk.ig <- read_excel(file.name.in2, sheet=2)  #PD data - Ig
  #datapk.tox <- read_excel(file.name.in2, sheet=3)  #toxicities
  #datapk.clin <- read_excel(file.name.in2, sheet=4)  #response
  #datapk.sero <- read_excel(file.name.in2, sheet=5)  #serotypes
  datapk.dose <- read_excel(file.name.in2, sheet=6)  #actual dose times

#-------------------------------------------------------------------------------
### Column names
##datanew
#pkd
  names(datanew.pkd) <- print(str_replace_all(names(datanew.pkd),"[ ()#]","."))
  str(datanew.pkd)
#act
  names(datanew.act) <- print(str_replace_all(names(datanew.act)," ","."))
  sort(names(datanew.act))  #two copies of each
  str(datanew.act)

##datapk
#demog
  names(datapk.demog) <- print(str_replace_all(names(datapk.demog)," ","."))
  str(datapk.demog)
#dose
  names(datapk.dose)
  str(datapk.dose)

#-------------------------------------------------------------------------------
#Plot PK data
#8 data points per ID, with 1 NA
  with(datanew.pkd, table(Patient.Number, useNA = "always"))

#2 lena dose levels, with 2 odd values
  with(datapk.dose, table(DOSE_LEVEL, useNA = "always"))
  with(datapk.dose, table(DRUG,DOSE_LEVEL))
  with(datapk.dose, table(SEQUENCE_NO_,DOSE_LEVEL))

#Number of patients (-1 for NAs)
  npat <- length(unique(datanew.pkd$Patient.Number))-1
  npat
  length(unique(datapk.dose$SEQUENCE_NO_))  #more patients in dose dataset than concentration dataset

#1 patient with concentrations that does not exist in dosage data
  unique(datanew.pkd$Patient.Number)[!unique(datanew.pkd$Patient.Number) %in% unique(datapk.dose$SEQUENCE_NO_)]

#60 mdv values
  with(datanew.pkd, table(Lenalidomide.Conc..nM., useNA = "always"))
  #all mdv values are on NA DV values
  temp <- with(datanew.pkd, table(Lenalidomide.Conc..nM., useNA = "always"))
  temp[2,dim(temp)[2]]




#-------------------------------------------------------------------------------
#Prepare continuous covariates
  names(demog.cont) <- c("ID","HT","WT","BSA","SECR","SECR")
  demog.cont$ID <- 1:length(demog.cont$ID)
  demog.cont2 <- rbind(data.frame(demog.cont[1:5],"DAY"=1), data.frame(demog.cont[c(1:4,6)],"DAY"=5))

  names(demog.cat)[c(1,9)] <- c("ID","RACE2")
  demog.cat$ID <- 1:length(demog.cat$ID)
  demog.cat2 <- demog.cat[c(1,9)]
  demog.cat2$RACE2[demog.cat$RACE2=="White"] <- 1
  demog.cat2$RACE2[demog.cat2$RACE2!=1] <- 2

#-------------------------------------------------------------------------------
#Convert datanew to old format

  datanew2 <- data.frame("ID" = datanew$ID, "STUDY" = 10016)

  #Cohort 1 -> 18 patients
							#INDUCTION - lenalidomide PO QD D1-D21 - cytarabine IV over 96hrs D5-D8 - idarubicin IV over 1hr D5-D7
						#CONSOLIDATION - lenalidomide PO QD D1-D14 - cytarabine continuous D5-D7 - idarubicin IV over 1hr D5-D6
  #Cohort 2	-> 14 patients
							#INDUCTION - lenalidomide PO QD D1-D21 - cytarabine IV over 24hrs D5-D11 - idarubicin IV over 1hr D5-D7
						#CONSOLIDATION - lenalidomide PO QD D1-D14 - cytarabine IV q24h D5,D7,D9

  datanew2$GRP <- datanew$Cohort

  #Dose Level
  #
  # -------------------   Cohort 1  	npat = 18
  #  1	25mg daily						npat = 12
  #  2  30mg daily		         		npat =  6
  # -------------------   Cohort 2	    npat = 14
  #  1  25mg daily    				   	npat =  7
  #  2  30mg daily      				npat =  7

  datanew2$DOSELVL <- 1
  datanew2$DOSELVL[datanew$dose..mg.==30] <- 2

  datanew2$DOSEMG <- datanew$dose..mg.

  datanew2$AMT <- datanew$amt				     #dose in mg

  datanew2$TIME <- datanew$time..hr.

  datanew2$DAY <- datanew$Day

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

  #datanew2$WT <- datanew$Weight..lbs./2.2		#conversion to kgs

  #datanew2$HT <- datanew$Height/3.28			#conversion to m

  #datanew2$BSA <- 0.007184*datanew2$WT**0.425*datanew2$HT**0.725

  #datanew2$BMI <- datanew2$WT/datanew2$HT**2



  #DXCATNUM == 2 -> Acute Myeloid Leukaemia		- Not Published NCT01132586
  datanew2$DXCATNUM <- 2

  #Caucasian 	1
  #??			2
  #??			3
  #with(datanew, table(Race, useNA = "always"))
  #datanew2$RACE <- datanew$Race

  #datanew2$RACE2 <- 2
  #datanew2$RACE2[datanew$Race == 1] <- 1
  #with(datanew2, table(RACE2, useNA = "always"))

  #datanew2$SECR <- datanew$SeCr..mg.dL.*88.4	#convert from mg/dL to umol/L

 #-----------------------------------------------------------------
	#Fixes for extra amt values
	temp <- rownames(datanew2[datanew2$TIME==0&!is.na(datanew2$DV)&!is.na(datanew2$AMT),])
	datanew2[temp,6] <- NA

  dataall <- merge(datanew2,demog.cont2)
  dataall <- merge(dataall,demog.cat2)[c("ID","STUDY","GRP","DOSELVL","DOSEMG","AMT","TIME","DAY","DV","MDV","LNDV","AGE","GEND","WT","HT","BSA","DXCATNUM","RACE2","SECR")]
  dataall <- orderBy(~ID+TIME, data=dataall)

#-------------------------------------------------------------------------------
# Check subject numbers
  with(dataall, table(ID))

  # Check the dose columns
  with(dataall, table(DOSEMG))
  with(dataall, table(DOSEMG,DOSELVL))	#dose by dose group

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
#Count missing covariate data
#Missing by Study

  covnames <- as.formula("~AGE+GEND+WT+HT+BSA+DXCATNUM+RACE2+SECR")
  covdata <- subset(dataall, select=c("ID","GRP","AGE","GEND","WT","HT","BSA","DXCATNUM","RACE2","SECR"))

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

  #with(dataallone,table(RATE,useNA="always"))

  #Sets missing to NA - use for continuous summary
  covdataone <- lapplyBy(~ID, data=covdata,  oneperID)
  covdataone <- bind.list(covdataone)
  dim(covdataone)


  #dataallone$RACEf <- factor(dataallone$RACE2, labels=c("White or Caucasian","Other"))

  dataallone$SEXf <- as.factor(dataallone$GEND)
  levels(dataallone$SEXf) <- c("female","male")

  dataallone$IDf <- as.factor(dataallone$ID)

  dataallone$DXCAT2f <- as.factor(dataallone$DXCATNUM)
  levels(dataallone$DXCAT2f) <- c("AML")

  dataallone$RACE2f <- as.factor(dataallone$RACE2)
  levels(dataallone$RACE2f) <- c("white", "non-white")

  #DXCAT2 "Acute Myeloid Leukaemia" <- 2

#-------------------------------------------------------------------------------
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
  #RATEtest <- summaryBy(RATE ~ ID, data=dataall, FUN=mean, na.rm=T)
  #RATEtest <- RATEtest[is.na(RATEtest$RATE.mean)==T,]
  #RATEtestID <- RATEtest$ID
  #RATEtestID
  #with(dataallone,table(RATE,useNA="always"))
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
  DVcountdata <-  ddply(dataall, .(GRP,ID,DAY), function(df) DV.present(df$DV))


  withDVbyGRPWEEK <- ddply(DVcountdata, .(GRP,DAY), function(df) sum(df$V1))  #GOLD
  withDVbyGRPWEEK


  #Not all subjects with Day 1 data also have Day 4 data
  filename.out <- paste(output.dir,"DVwith_group_week.csv",sep="/")
  write.csv(withDVbyGRPWEEK, file=filename.out)

#-------------------------------------------------------------------------------
#Subset some plot data

  plotdata <- subset(dataall)
  BINnumber <- 3

  plotdata$DOSEMGf <- as.factor(plotdata$DOSEMG)

  plotdata$GRPf <- as.factor(plotdata$GRP)
  levels(plotdata$STUDYf) <- paste("Group",levels(plotdata$STUDYf))

  #plotdata$VOSf <- as.factor(plotdata$VOS)
  #levels(plotdata$VOSf) <- c("Placebo","Adrug")

  plotdata$DAYf <- as.factor(plotdata$DAY)
  levels(plotdata$DAYf) <- paste("Day",levels(plotdata$DAYf))

  plotdata$SEXf <- as.factor(plotdata$GEND)
  levels(plotdata$SEXf) <- c("female","male")

  #plotdata$RACEf <- as.factor(plotdata$RACE2)
  #levels(plotdata$RACEf) <- c("White or Caucasian","Other")

  plotdata$DOSELVLf <- as.factor(plotdata$DOSELVL)
  levels(plotdata$DOSELVLf) <- paste("Dose Level",levels(plotdata$DOSELVLf))

  #plotdata$DOSE_bin <- cut2(plotdata$DOSEMG, g=BINnumber)

  plotdata$AGE_bin <- cut2(plotdata$AGE, g=BINnumber)

  #plotdata$WT_bin <- cut2(plotdata$WT, g=BINnumber)

  #plotdata$HT_bin <- cut2(plotdata$HT, g=BINnumber)

  #plotdata$BSA_bin <- cut2(plotdata$BSA, g=BINnumber)

  #plotdata$BMI_bin <- cut2(plotdata$BMI, g=BINnumber)

  plotdata$DXCAT2f <- as.factor(plotdata$DXCATNUM)
  levels(plotdata$DXCAT2f) <- c("Acute Myeloid Leukaemia")

  filename.out <- paste(output.dir,"plotdata.csv",sep="/")
  write.csv(plotdata, file=filename.out, row.names=FALSE)

#-------------------------------------------------------------------------------
#Basic PK plots

 #Conc vs Time
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=plotdata)  #, colour=AMTMGf
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DAYf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)")  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Day")
  plotobj

  filename.out <- paste(output.dir,"ConcObs_vs_TIME_by_DAY",sep="/")
  to.png(plotobj,filename.out)

 #Conc vs TIME
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~DAY)
  plotobj

  filename.out <- paste(output.dir,"ConcObs_vs_TIME_by_DOSELVL",sep="/")
  to.png(plotobj,filename.out)

 #Conc vs TIME per ID
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DAYf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Day")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"ConcObs_vs_TIME_by_ID",sep="/")
  to.png(plotobj,filename.out)

#DVNORM vs TIME
  plotobj <- NULL
  titletext <- paste("Dose Normalised Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DVNORM, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~DAY)
  plotobj

  filename.out <- paste(output.dir,"DVNORM_vs_TIME",sep="/")
  to.png(plotobj,filename.out)


#-------------------------------------------------------------------------------
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
  plotByFactor("DAYf","Day")
  #plotByFactor("DOSE_bin","Binned Dose (mg)")
  plotByFactor("AGE_bin","Binned Age (years)")
  #plotByFactor("WT_bin","Binned Weight (kg)")
  #plotByFactor("HT_bin","Binned Height (kg)")
  #plotByFactor("BSA_bin","Binned BSA (m2)")
  #plotByFactor("BMI_bin","Binned BMI (kg per m2)")
  #plotByFactor("DXCAT2f","Disease category")


#-------------------------------------------------------------------------------
#Summarize Categorical covariates

  dataallone$SUMCOL <- "All Groups"

  #Sex
  #covCatSexStudy <- with(dataallone,ftable(GRP,SEXf, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatSexStudy  <- data.frame("COV"="SEX",covCatSexStudy)


  #Race
  #covCatRaceStudy <- with(dataallone,ftable(GRP,RACEf, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatRaceStudy  <- data.frame("COV"="RACE",covCatRaceStudy)


  #DXcategory
  #covCatDXcatStudy <- with(dataallone,ftable(GRP,DXCAT2f, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatDXcatStudy <- data.frame("COV"="DXCAT",covCatDXcatStudy)


  #Collate
  #covCatTable <- rbind(covCatSexStudy,covCatRaceStudy,covCatDXcatStudy)
  #covCatTable

  #Return to original order
  #covCatTable <- orderBy(~COV, covCatTable)  #GOLD - sort by factor levels

  #Define reassignment of column names (if any) rtf allowed
  #colNames <- c(
  #              "COV","Covariate Code",
  #              "CATEGORY","Category",
  #              "RACE","Race",
  #              "SEX","Gender",
  #              "Freq","Count",
	#			"DXCAT","Disease"
                )

  #colData <- data.frame(matrix(colNames, byrow=T, ncol=2),stringsAsFactors=F)
  #names(colData) <- c("ColNameOld","ColNameNew")

  #Reassign column names
  #covCatTable$COV <- gsub.all(colData$ColNameOld,colData$ColNameNew,covCatTable$COV)
  #names(covCatTable) <- gsub.all(colData$ColNameOld,colData$ColNameNew,names(covCatTable))

  #Write as an rtf file
  #filename.out <- paste(output.dir,"CatCovSummary",sep="/")
  #write.csv(covCatTable, file=filename.out)

#-------------------------------------------------------------------------------
  #Pairs-plot of continuous covariates
  plotobj <- NULL
  covdataonecont <- subset(covdataone, select=c("AGE","WT","BSA","BMI","SECR"))
  plotobj <- ggpairs(na.omit(covdataonecont))

  filename.out <- paste(output.dir,"Overview_cont_cov_pairs",sep="/")
  to.png(plotobj,filename.out)


#-------------------------------------------------------------------------------
 #Pairs-plot of categorical covariates

  covcatdataonecat <- subset(dataallone, select=c("SEXf","RACEf","DXCAT2f"))

  plotobj <- NULL
  plotobj <- ggpairs(na.omit(covcatdataonecat))

  filename.out <- paste(output.dir,"Overview_cat_cat_pairs",sep="/")
  to.png(plotobj,filename.out)


#-------------------------------------------------------------------------------
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
  covdataone$AGEBIN <- cut2(covdataone$AGE, cuts=c(18,50,75,85))
  #covdataone$AGEBIN <- as.numeric(paste(covdataone$AGEBIN))


  #GEND Summary
  #with(covdataone, table(GEND))

  #with(covdataone, table(GEND,AGEBIN))


  #RACE Summary
  #RACEtable <-  with(dataallone, table(RACEf))

  #filename.out <- paste(output.dir,"RACEtable.csv",sep="/")
  #write.csv(RACEtable, file=filename.out, row.names=T)


#-------------------------------------------------------------------------------
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
  #plotIndexCont("WT","Weight~(kg)")
  #plotIndexCont("BSA","Body~Surface~Area~(m^2)")
  #plotIndexCont("BMI","Body~Mass~Index~(kg/m^2)")
  plotIndexCat("GEND","Patient~Sex")
  #plotIndexCat("RACE2","Patient~Race")
  #plotIndexCat("DXCATNUM","Diagnosis~Category")
  #plotIndexCont("SECR","Serum~Creatinine~(umol/L)")

	#Data prep
	# [1] "#ID"      "STUDY"    "XSAMP"    "GRP"      "DOSELVL"  "DOSEMG"   "AMT"      "RATE"     "TIME"
	#[10] "TAD"      "DAY"      "DV"       "MDV"      "LNDV"     "AGE"      "GEND"     "WT"       "HT"
	#[19] "BSA"      "BMI"      "DXCATNUM" "RACE"     "RACE2"    "SECR"     "DVNORM"   "ADDL"     "II"

	dataFIX <- data.frame("ID" = (dataall$ID+121), "STUDY" = dataall$STUDY)

  dataFIX$XSAMP <- 0
  dataFIX$GRP <- dataall$GRP+6

  dataFIX$DOSELVL <- dataall$DOSELVL

  dataFIX$DOSEMG <- dataall$DOSEMG
	dataFIX$AMT <- dataall$AMT
  dataFIX$RATE <- 0
  dataFIX$TIME <- dataall$TIME+(dataall$DAY-1)*24
	dataFIX$TAD <- dataall$TIME

  dataFIX$DAY <- dataall$DAY
  dataFIX$DV <- dataall$DV
  dataFIX$MDV <- dataall$MDV
  dataFIX$LNDV <- dataall$LNDV
  dataFIX$AGE <- dataall$AGE
  dataFIX$GEND <- dataall$GEND

  dataFIX$WT <- dataall$WT
  dataFIX$HT <- dataall$HT
  dataFIX$BSA <- dataall$BSA
  dataFIX$BMI <- dataall$WT/dataall$HT**2
  dataFIX$DXCATNUM <- dataall$DXCATNUM
  dataFIX$RACE <- dataall$RACE2
  dataFIX$RACE2 <- dataall$RACE2
  dataFIX$SECR <- dataall$SECR*88.4	#convert from mg/dL to umol/L
  dataFIX$DVNORM <- dataall$DVNORM

  dataFIX$ADDL <- NA
  dataFIX$ADDL[!is.na(dataFIX$AMT)] <- 20

  dataFIX$II <- NA
	dataFIX$II[!is.na(dataFIX$AMT)] <- 24

	dataFIX[(dataFIX$DAY==5&!is.na(dataFIX$AMT)),] <- NA
	dataFIX <- orderBy(~ID+TIME+DAY+AMT, data=dataFIX)
	colnames(dataFIX)[1] <- "#ID"

	filename.out <- paste(output.dir,"10016_finaldata.csv",sep="/")
  write.csv(dataFIX[1:(dim(dataFIX)[1]-25),], file=filename.out, row.names=FALSE)


#------------------
#Covariate data
	# [1] "UID"      "ID"       "STUDY"    "GRP"      "DOSELVL"  "DOSEMG"   "AGE"      "GEND"
	# [8] "WEIGHTLB" "HEIGHTFT" "DXCATNUM" "RACE"     "SECRMGDL"

	dataCOV <- data.frame("UID" = (dataallone$ID+121), "ID" = dataallone$ID, "STUDY" = dataallone$STUDY)

  dataCOV$GRP <- dataallone$GRP+6
  dataCOV$DOSELVL <- dataallone$DOSELVL
  dataCOV$DOSEMG <- dataallone$DOSEMG
  dataCOV$AGE <- dataallone$AGE
  dataCOV$GEND <- dataallone$GEND
  dataCOV$WEIGHTLB <- dataallone$WT*2.2
  dataCOV$HEIGHTFT <- dataallone$HT*3.28
  dataCOV$DXCATNUM <- dataallone$DXCATNUM
  dataCOV$RACE <- dataallone$RACE2
  dataCOV$SECRMGDL <- dataallone$SECR

	filename.out <- paste(output.dir,"10016_covdata.csv",sep="/")
  write.csv(dataCOV, file=filename.out, row.names=FALSE)
