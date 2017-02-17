###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_clin_05115"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)
  library(GGally)
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
  file.name.in <- "RAW_Clinical/rawdata-lena_05115_allinfo_creatinine.csv"
  datanew <- read.csv(file.name.in, stringsAsFactors=F, na.strings=c("."))

  ## Wt_Ht_BMI.xlsx - 13 sheets of data
  file.name.in3 <- "RAW_Clinical/rawdata-lena_05115_06003_Wt_Ht_BMI.xlsx"
  demog.ht <- read_excel(file.name.in3, sheet=1)  #demographics

  ##corr_analysis.xls - 9 sheets of data
  file.name.in2 <- "RAW_Clinical/rawdata-lena_05115_corr_analysis.xls"
  datapd.in <- data.frame(read_excel(file.name.in2, sheet=9)[-1,-(2:24)])  #rawdata
  names(datapd.in) <- datapd.in[1,]
  dataraceone <- data.frame(
    "ID" = as.numeric(datapd.in[2:21,1]),
    "RACE" = ifelse(datapd.in[2:21,6] == "B", 2, 1))

#------------------------------------------------------------------------------------
#Column names
  #As presented
  names(datanew)

  #Sorted
  sort(names(datanew))

  #Structure
  str(datanew)

#------------------------------------------------------------------------------------
#Plot PK data

#9 data entries for each individual
  with(datanew, table(ID, useNA = "always"))

#3 dose levels
  with(datanew, table(dose, useNA = "always"))
  with(datanew, table(dose,ID))

#23 mdv values
  with(datanew, table(mdv, useNA = "always"))
  #all mdv values are on NA DV values
  temp <- with(datanew, table(mdv, dv..ug.mL., useNA = "always"))
  temp[2,dim(temp)[2]]

  #Number of patients
  npat <- length(unique(datanew$ID))
  npat

#-------------------------------------------------------------------------------
#Get heights from demog.ht
#Column names
  #As presented
  names(demog.ht)

  data.ht <- demog.ht[1:21,c(1,3,11)]
  names(data.ht) <- c("ID","RACE","HT")
  data.ht$ID <- 1:21
  data.ht <- data.ht[!is.na(data.ht$HT),]
  data.ht$RACE <- ifelse(data.ht$RACE=="W", 1, 2)
  datanew <- merge(datanew,data.ht)

#-------------------------------------------------------------------------------
#Convert datanew to old format

  datanew2 <- data.frame("ID" = datanew$ID, "STUDY" = 5115)

  #Dose Level 											     		   total patients = 20
  #
  # -----------------   npat = 20
  #	 1	15mg daily			npat =  1
  #  2	20mg daily			npat =  6
  #  3  25mg daily      npat = 13

  datanew2$DOSELVL <- 3
  datanew2$DOSELVL[datanew$dose==20] <- 2
  datanew2$DOSELVL[datanew$dose==15] <- 1

  datanew2$DOSEMG <- datanew$dose

  datanew2$AMT <- datanew$amt				     #dose in mg

  datanew2$TIME <- datanew$time..hr.

  datanew2$EVID <- datanew$dvid

  datanew2$DV <- datanew$dv..ug.mL.    			 #ug/ml
  datanew2$DV[datanew2$DV <= 0] <- NA

  datanew2$MDV <- datanew$mdv
  datanew2$MDV[datanew2$DV<0.00025926] <- 1

  #datanew2$BLQ <- 0
  #datanew2$BLQ[datanew$NOTE == "DV_BLQ"] <- 1
  #with(datanew2, table(BLQ, useNA = "always"))
  #datanew2$DV[datanew2$BLQ==1] <- NA

  datanew2$LNDV <- log(datanew2$DV)

  #datanew2$DNUM <- datanew$Dose.no
  #datanew2$DNUM <- unlist(lapplyBy(~ID, data=datanew2, function(d) impute(d$DNUM)))

  #datanew2$OCC <- datanew2$DNUM

  datanew2$AGE <- datanew$age

  datanew2$GEND <- datanew$sex				  	#1 is male, 0 is female
  with(datanew2, table(GEND, useNA = "always"))

  datanew2$WT <- datanew$wt..lbs./2.2		#conversion to kgs

  datanew2$HT <- as.numeric(datanew$HT)*100

  #datanew2$BSA <- 0.007184*datanew2$WT**0.425*datanew2$HT**0.725

  #datanew2$BMI <- datanew2$WT/datanew2$HT**2

  #DXCATNUM == 4 -> Relapsed Multiple Myeloma - Hofmeister et. al 2011
  datanew2$DXCATNUM <- 4


  #Caucasian     	1
  #Non-Caucasian 	2
  datarace <- merge(datanew2, dataraceone)
  datanew2$RACE <- datarace$RACE
  datanew2$RACE2 <- datarace$RACE

  datanew2$SECR <- datanew$Cr..mg.dL.*88.4	#convert from mg/dL to umol/L

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

  covnames <- as.formula("~AGE+GEND+WT+DXCATNUM+SECR")
  covdata <- subset(dataall, select=c("ID","DOSELVL","AGE","GEND","WT","HT","DXCATNUM","SECR"))

  #Reassign missing
  covdata[covdata==-1] <- NA

  #finish off
  missingbystudy <- ddply(covdata, .(DOSELVL), colwise(calculate.percent.missing))

  filename.out <- paste(output.dir,"Missing_by_Group.csv",sep="/")
  write.csv(missingbystudy, file=filename.out, row.names=F)



#Missing by Subject
  missingbysubject <- ddply(covdata, .(DOSELVL,ID), colwise(calculate.percent.missing))
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
  levels(dataallone$DXCAT2f) <- c("CLL")

  #DXCAT2 "Relapsed Multiple Myeloma" <- 4

#-----------------------------------------------------------------------
#Summary of study characteristics

 #Do all subjects have PK data
  DVtest <- summaryBy(DXCATNUM ~ ID, data=dfnew, FUN=mean, na.rm=T)
  DVtest <- DVtest[is.na(DVtest$DXCATNUM.mean)==T,]
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


 #DV count by Study
 #Calculates data for Report Table 1
  DVcount <- summaryBy(DV ~ DOSELVL, data=dataall, FUN=lengthNA)
  names(DVcount) <- c("Group","DVcount")
  DVcount


 #Subject count by Study
  SUBcount <- ddply(dataall, .(DOSELVL), function(df) count.unique(df$ID))
  names(SUBcount) <- c("Group","SUBcount")
  SUBcount

 #Dose count by Study
  AMTcount <- ddply(dataall, .(DOSELVL), function(df) lengthNA(df$AMT))
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
  #DVcount <- summaryBy(DV ~ GRP+DOSELVL, data=dataall, FUN=lengthNA)
  #names(DVcount) <- c("Group","Dose Level","DVcount")
  #DVcount


#Subject count by Group and DoseLevel
  #SUBcount <- ddply(dataall, .(GRP,DOSELVL), function(df) count.unique(df$ID))
  #names(SUBcount) <- c("Group","Dose Level","SUBcount")
  #SUBcount


#Dose count by Group and DoseLevel
  #AMTcount <- ddply(dataall, .(GRP,DOSELVL), function(df) lengthNA(df$AMT))
  #names(AMTcount) <- c("Group","Dose Level","AMTcount")
  #AMTcount


#Average DV and AMT per Subject
  #DVsum <- cbind(DVcount,SUBcount,AMTcount)
  #DVsum$DVperSUB <-  round(DVsum$DVcount/DVsum$SUBcount,0)
  #DVsum$AMTperSUB <-  round(DVsum$AMTcount/DVsum$SUBcount,0)
  #DVsum

  #filename.out <- paste(output.dir,"DVsum_DOSELVL.csv",sep="/")
  #write.csv(DVsum, file=filename.out)

#DV data present by Group
  DV.present <- function(x)  if (any(is.numeric(x))==T) result <- 1 else result <- 0
  DVcountdata <-  ddply(dataall, .(DOSELVL,ID), function(df) DV.present(df$DV))


  withDVbyGRP <- ddply(DVcountdata, .(DOSELVL), function(df) sum(df$V1))  #GOLD
  withDVbyGRP


  #Not all subjects with Day 1 data also have Day 4 data
  filename.out <- paste(output.dir,"DVwith_group.csv",sep="/")
  write.csv(withDVbyGRP, file=filename.out)

#----------------------------------------------------------------------------------------------------------------------
#Subset some plot data

  plotdata <- subset(dataall)
  BINnumber <- 3

  plotdata$DOSEMGf <- as.factor(plotdata$DOSEMG)

  #plotdata$GRPf <- as.factor(plotdata$GRP)
  #levels(plotdata$STUDYf) <- paste("Group",levels(plotdata$STUDYf))

  #plotdata$VOSf <- as.factor(plotdata$VOS)
  #levels(plotdata$VOSf) <- c("Placebo","Adrug")

  #plotdata$WEEKf <- as.factor(plotdata$WEEK)
  #levels(plotdata$WEEKf) <- paste("Week",levels(plotdata$WEEKf))

  plotdata$SEXf <- as.factor(plotdata$GEND)
  levels(plotdata$SEXf) <- c("female","male")

  #plotdata$RACEf <- as.factor(plotdata$RACE2)
  #levels(plotdata$RACEf) <- c("White or Caucasian","Other")

  plotdata$DOSELVLf <- as.factor(plotdata$DOSELVL)
  levels(plotdata$DOSELVLf) <- paste("Dose Level",levels(plotdata$DOSELVLf))

  #plotdata$DOSE_bin <- cut2(plotdata$DOSEMG, g=BINnumber)

  plotdata$AGE_bin <- cut2(plotdata$AGE, g=BINnumber)

  plotdata$WT_bin <- cut2(plotdata$WT, g=BINnumber)

  plotdata$HT_bin <- cut2(plotdata$HT, g=BINnumber)

  #plotdata$BSA_bin <- cut2(plotdata$BSA, g=BINnumber)

  #plotdata$BMI_bin <- cut2(plotdata$BMI, g=BINnumber)

  plotdata$DXCAT2f <- as.factor(plotdata$DXCATNUM)
  levels(plotdata$DXCAT2f) <- c("Relapsed Multiple Myeloma")

  filename.out <- paste(output.dir,"plotdata.csv",sep="/")
  write.csv(plotdata, file=filename.out, row.names=FALSE)

#----------------------------------------------------------------------------------------------------------------------
#Basic PK plots

 #Conc vs TIME Week 1
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj

  filename.out <- paste(output.dir,"ConcObs_vs_TIME_by_DOSELVL",sep="/")
  to.png(plotobj,filename.out)

 #Conc vs TIME Week 1 per ID
  plotobj <- NULL
  titletext <- paste("Observed Concentrations\n")
  plotobj <- ggplot(data=subset(plotdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
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
  plotobj

  filename.out <- paste(output.dir,"DVNORM_vs_TIME",sep="/")
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
  #plotByFactor("GRPf","Group")
  plotByFactor("SEXf","SEX")
  plotByFactor("DOSELVLf","Dose Level")
  #plotByFactor("WEEKf","Week")
  #plotByFactor("DOSE_bin","Binned Dose (mg)")
  plotByFactor("AGE_bin","Binned Age (years)")
  plotByFactor("WT_bin","Binned Weight (kg)")
  plotByFactor("HT_bin","Binned Height (kg)")
  #plotByFactor("BSA_bin","Binned BSA (m2)")
  #plotByFactor("BMI_bin","Binned BMI (kg per m2)")
  #plotByFactor("DXCAT2f","Disease category")


#-------------------------------------------------------------------------
#Summarize Categorical covariates

  #dataallone$SUMCOL <- "All Groups"

  #Sex
  #covCatSexStudy <- with(dataallone,ftable(DOSELVL,SEXf, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatSexStudy  <- data.frame("COV"="SEX",covCatSexStudy)


  #Race
  #covCatRaceStudy <- with(dataallone,ftable(DOSELVL,RACEf, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatRaceStudy  <- data.frame("COV"="RACE",covCatRaceStudy)


  #DXcategory
  #covCatDXcatStudy <- with(dataallone,ftable(DOSELVL,DXCAT2f, useNA="ifany", dnn=c("GRP","CATEGORY")))
  #covCatDXcatStudy <- data.frame("COV"="DXCAT",covCatDXcatStudy)


  #Collate
  #covCatTable <- rbind(covCatSexStudy,covCatRaceStudy,covCatDXcatStudy)
  #covCatTable

  #Return to original order
  #covCatTable <- orderBy(~COV, covCatTable)  #GOLD - sort by factor levels

  #Define reassignment of column names (if any) rtf allowed
  #colNames <- c(
  #              "COV","Covariate Code",
#                "CATEGORY","Category",
#                "RACE","Race",
#                "SEX","Gender",
#                "Freq","Count",
#				"DXCAT","Disease"
#                )

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
  plotobj <- NULL
  covdataonecont <- subset(covdataone, select=c("AGE","WT","SECR"))
  plotobj <- ggpairs(na.omit(covdataonecont))

  filename.out <- paste(output.dir,"Overview_cont_cov_pairs",sep="/")
  to.png(plotobj,filename.out)

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
  #plotIndexCont("BSA","Body~Surface~Area~(m^2)")
  #plotIndexCont("BMI","Body~Mass~Index~(kg/m^2)")
  plotIndexCat("GEND","Patient~Sex")
  #plotIndexCat("RACE2","Patient~Race")
  #plotIndexCat("DXCATNUM","Diagnosis~Category")
  plotIndexCont("SECR","Serum~Creatinine~(umol/L)")

#--------------------
#Data prep
	# [1] "#ID"      "STUDY"    "XSAMP"    "GRP"      "DOSELVL"  "DOSEMG"   "AMT"      "RATE"     "TIME"
	#[10] "TAD"			 "DAY"      "DV"       "MDV"      "LNDV"     "AGE"      "GEND"     "WT"       "HT"
	#[19] "BSA"      "BMI"      "DXCATNUM" "RACE"     "RACE2"    "SECR"     "DVNORM"   "ADDL"     "II"

  dataFIX <- data.frame("ID" = (dataall$ID+75), "STUDY" = dataall$STUDY)

  dataFIX$XSAMP <- 0
  dataFIX$GRP <- 4
  dataFIX$DOSELVL <- dataall$DOSELVL
  dataFIX$DOSEMG <- dataall$DOSEMG
	dataFIX$AMT <- dataall$AMT
  dataFIX$RATE <- 0
  dataFIX$TIME <- dataall$TIME
	dataFIX$TAD <- dataall$TIME
  dataFIX$DAY <- 1
  dataFIX$DV <- dataall$DV
  dataFIX$MDV <- dataall$MDV
  dataFIX$LNDV <- dataall$LNDV
  dataFIX$AGE <- dataall$AGE
  dataFIX$GEND <- dataall$GEND
  dataFIX$WT <- dataall$WT

  dataFIX$HT <- dataall$HT
  dataFIX$BSA <- NA
  dataFIX$BMI <- NA
  dataFIX$DXCATNUM <- dataall$DXCATNUM
  dataFIX$RACE <- dataall$RACE
  dataFIX$RACE2 <- dataall$RACE
  dataFIX$SECR <- dataall$SECR
  dataFIX$DVNORM <- dataall$DVNORM

  dataFIX$ADDL <- NA
  dataFIX$ADDL[!is.na(dataFIX$AMT)] <- 20

  dataFIX$II <- NA
	dataFIX$II[!is.na(dataFIX$AMT)] <- 24

	colnames(dataFIX)[1] <- "#ID"

	filename.out <- paste(output.dir,"05115_finaldata.csv",sep="/")
  write.csv(dataFIX, file=filename.out, row.names=FALSE)

#------------------
#Covariate data
	# [1] "UID"      "ID"       "STUDY"    "GRP"      "DOSELVL"  "DOSEMG"   "AGE"      "GEND"
	# [8] "WEIGHTLB" "HEIGHTFT" "DXCATNUM" "RACE"     "SECRMGDL"

	dataCOV <- data.frame("UID" = (dataallone$ID+75), "ID" = dataallone$ID, "STUDY" = dataallone$STUDY)

  dataCOV$GRP <- 4
  dataCOV$DOSELVL <- dataallone$DOSELVL
  dataCOV$DOSEMG <- dataallone$DOSEMG
  dataCOV$AGE <- dataallone$AGE
  dataCOV$GEND <- dataallone$GEND
  dataCOV$WT <- dataallone$WT
  dataCOV$HT <- dataallone$HT
  dataCOV$DXCATNUM <- dataallone$DXCATNUM
  dataCOV$RACE <- dataallone$RACE
  dataCOV$SECRMGDL <- dataallone$SECR/88.4

	filename.out <- paste(output.dir,"05115_covdata.csv",sep="/")
  write.csv(dataCOV, file=filename.out, row.names=FALSE)

# ------------------------------------------------------------------------------
# Clean Data (non-nmprep)
  library(dplyr)
  clean.data <- dataFIX
  names(clean.data)[1] <- "ID"
  clean.data <- clean.data %>%
    select(c(ID, STUDY, DOSEMG, TIME, TAD, DAY, DV, MDV,
      AGE, GEND, WT, HT, DXCATNUM, RACE, SECR)) %>%
    rename(DXCAT = DXCATNUM, SEX = GEND, DVMGL = DV, WTKG = WT, HTCM = HT, SECRUMOLL = SECR)
  clean.data$ID <- clean.data$ID - 75

  dxcat.l <- "MM"
  clean.data$DXCAT <- factor(clean.data$DXCAT, levels = c(4))
  levels(clean.data$DXCAT) <- dxcat.l

  sex.l <- c("F", "M")
  clean.data$SEX <- factor(clean.data$SEX, levels = c(0, 1))
  levels(clean.data$SEX) <- sex.l

  race.l <- c("W", "B")
  clean.data$RACE <- factor(clean.data$RACE, levels = c(1, 2))
  levels(clean.data$RACE) <- race.l

  filename.out <- paste(output.dir,"05115_cleandata.csv",sep="/")
  write.csv(dataCOV, file=filename.out, row.names=FALSE)
