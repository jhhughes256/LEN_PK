#SCRIPT TO GENERATE A VPC
#New Style VPCs

#Before you start need to do the following steps:
	#1 - Type into NONMEM "NMCTL (Model file name.ctl)"
	     #Only do this if you want the final parameter estimates to replace the initial estimates;
		 #Otherwise duplicate the ctl and rename the copy as  Filename_VPC.ctl - Continue now from step 4
	#2 - NONMEM generates a copy of   this ctl with the same name but final parameter estimates added and renames the orginal .ctl as a .org
	#3 - Rename the newly generated .ctl as filename_VCP.ctl and change the orginal (.org) back to its ctl form
	#4 - In the Filename_VPC.ctl change it to be $SIMULATION (1234567)  ONLYSIM SUBPROBLEMS=500
	#5 - Comment out $ESTIMATION and $COVARIANCE
	#5 - Run this Filename_VPC.ctl and process the results/fit file with this r script

#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()


#Set the working directory - to the parent directory where you do the modelling/VPC
	master.dir<-"E:/Hughes/Data/PK/"
	setwd(master.dir)

#Load libraries
	library(R2HTML)
	library(ggplot2)
	library(doBy)
	library(stringr)
	library(Hmisc)
	library(grid)
	library(plyr)
	library(reshape2)
	library(npde)

#Source WfNviaR functions file - Need to have this loaded as it runs most of the processing
	source("E:/Hughes/functions_utility.r")

#Use custom ggplot 2 theme
	theme_custom <- theme_set(theme_bw(18))
	theme_custom <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"))

#-------------------------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
theme_bw2 <- theme_set(theme_bw(base_size = 22))
theme_bw2 <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"),
                         axis.title.x=element_text(size = 18, vjust = 0),
                         axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
                         strip.text.x=element_text(size = 16),
                         strip.text.y=element_text(size = 16, angle = 90))

# ------------------------------------------------------------------------------
	runname <- c(
		"RUN058_DES_1C8TAwAP2_PPV_CORCLVKA_FFM_VPC",
		"RUN025_EXTVAL10_RUV1_ALLOM_VPC",
		"RUN029_EXTVAL_LOPEZ_VPC")

	simdata1 <- read.csv(paste0("REDO/", runname[1], ".nm7/", runname[1], ".fit.csv"), stringsAsFactors=F, na.strings=".")
	simdata1 <- subset(simdata1, MDV == 0)
  simdata2 <- read.csv(paste0("NO_ALLOM/", runname[2], ".nm7/", runname[2], ".fit.csv"), stringsAsFactors=F, na.strings=".")
	simdata2 <- subset(simdata2, MDV == 0)
  simdata3 <- read.csv(paste0("NO_ALLOM/", runname[3], ".nm7/", runname[3], ".fit.csv"), stringsAsFactors=F, na.strings=".")
	simdata3 <- subset(simdata3, MDV == 0)

#-------------------------------------------------------------------------------
#Read the original data
	orgdata1 <- read.csv("REDO/nmprep_allstudies.csv", stringsAsFactors=F, na.strings=".")
	orgdata1 <- rename(orgdata1, c("X.ID" = "ID"))
	orgdata1 <- subset(orgdata1, MDV == 0 & EVID <= 1)

	orgdata2 <- read.csv("NO_ALLOM/extval_10156.csv", stringsAsFactors=F, na.strings=".")
	orgdata2 <- rename(orgdata2, c("X.ID" = "ID"))
	orgdata2 <- subset(orgdata2, MDV == 0 & EVID <= 1)

	orgdata3 <- read.csv("NO_ALLOM/extval_allstudies.csv", stringsAsFactors=F, na.strings=".")
	orgdata3 <- rename(orgdata3, c("X.ID" = "ID"))
	orgdata3 <- subset(orgdata3, MDV == 0 & EVID <= 1)

 #------------------------------------------------------------------------------
#Assign some factors - ORG.data
	#Bin time after first dose
	orgdata1$TADBIN <- cut2(orgdata1$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02,22.02), levels.mean=T)
	orgdata1$TADBIN <- as.numeric(paste(orgdata1$TADBIN))

	orgdata2$TADBIN <- cut2(orgdata2$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02), levels.mean=T)
	orgdata2$TADBIN <- as.numeric(paste(orgdata2$TADBIN))

	orgdata3$TADBIN <- cut2(orgdata3$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02,22.02), levels.mean=T)
	orgdata3$TADBIN <- as.numeric(paste(orgdata3$TADBIN))

  # Assign some factors
  orgdata1$IDf <- as.factor(orgdata1$ID)
  orgdata1$RACEf <- factor(orgdata1$RACE, labels=c("White", "Other"))
  orgdata1$GENDf <- factor(orgdata1$GEND, labels=c("F", "M"))
  orgdata1$DXCATf <- factor(orgdata1$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
  orgdata1$STUDYf <- factor(orgdata1$STUDY)
	orgdata1$DOSEf <- factor(ifelse(orgdata1$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	orgdata1$CRCLf <- factor(ifelse(orgdata1$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	orgdata2$IDf <- as.factor(orgdata2$ID)
  orgdata2$RACEf <- factor(orgdata2$RACE, labels=c("White", "Other"))
  orgdata2$GENDf <- factor(orgdata2$GEND, labels=c("F", "M"))
  orgdata2$DXCATf <- factor(orgdata2$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
  orgdata2$STUDYf <- factor(orgdata2$STUDY)
	orgdata2$DOSEf <- factor(ifelse(orgdata2$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	orgdata2$CRCLf <- factor(ifelse(orgdata2$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	orgdata3$IDf <- as.factor(orgdata3$ID)
	orgdata3$RACEf <- factor(orgdata3$RACE, labels=c("White", "Other"))
	orgdata3$GENDf <- factor(orgdata3$GEND, labels=c("F", "M"))
	orgdata3$DXCATf <- factor(orgdata3$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
	orgdata3$STUDYf <- factor(orgdata3$STUDY)
	orgdata3$DOSEf <- factor(ifelse(orgdata3$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	orgdata3$CRCLf <- factor(ifelse(orgdata3$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

#-------------------------------------------------------------------------------
#Assign some factors - SIM.data
	#Bin time after first dose
	simdata1$TADBIN <- cut2(simdata1$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02,22.02), levels.mean=T)
	simdata1$TADBIN <- as.numeric(paste(simdata1$TADBIN))

	simdata2$TADBIN <- cut2(simdata2$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02), levels.mean=T)
	simdata2$TADBIN <- as.numeric(paste(simdata2$TADBIN))

	simdata3$TADBIN <- cut2(simdata3$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02,22.02), levels.mean=T)
	simdata3$TADBIN <- as.numeric(paste(simdata3$TADBIN))

	# Assign some factors
	simdata1$IDf <- as.factor(simdata1$ID)
	simdata1$RACEf <- factor(simdata1$RACE, labels=c("White", "Other"))
	simdata1$GENDf <- factor(simdata1$GEND, labels=c("F", "M"))
	simdata1$DXCATf <- factor(simdata1$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
	simdata1$STUDYf <- factor(simdata1$STUDY)
	simdata1$DOSEf <- factor(ifelse(simdata1$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	simdata1$CRCLf <- factor(ifelse(simdata1$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	simdata2$IDf <- as.factor(simdata2$ID)
	simdata2$RACEf <- factor(simdata2$RACE, labels=c("White", "Other"))
	simdata2$GENDf <- factor(simdata2$GEND, labels=c("F", "M"))
	simdata2$DXCATf <- factor(simdata2$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
	simdata2$STUDYf <- factor(simdata2$STUDY)
	simdata2$DOSEf <- factor(ifelse(simdata2$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	simdata2$CRCLf <- factor(ifelse(simdata2$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	simdata3$IDf <- as.factor(simdata3$ID)
	simdata3$RACEf <- factor(simdata3$RACE, labels=c("White", "Other"))
	simdata3$GENDf <- factor(simdata3$GEND, labels=c("F", "M"))
	simdata3$DXCATf <- factor(simdata3$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
	simdata3$STUDYf <- factor(simdata3$STUDY)
	simdata3$DOSEf <- factor(ifelse(simdata3$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	simdata3$CRCLf <- factor(ifelse(simdata3$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

#-------------------------------------------------------------------------------
#Confidence intervals - from function utility

	CI90lo <- function(x) quantile(x, probs=0.05)
	CI90hi <- function(x) quantile(x, probs=0.95)

	CI95lo <- function(x) quantile(x, probs=0.025)
	CI95hi <- function(x) quantile(x, probs=0.975)

	#-----------------------------------------------------------------------------

	#Change working directory, so that plots are stored within the VPC folder
	setwd("E:/Hughes/Data/PK/Plots")

#-------------------------------------------------------------------------------
###GENERATE A pcVPC
###Bergstrand et al 2011 - Prediction-Corrected Visual Predictive Checks for Diagnosing Nonlinear Mixed-Effects Models
# ------------------------------------------------------------------------------
# Simdata1

#Calculate the median PRED for each TADBIN
	simdata1$PRED <- as.numeric(simdata1$PRED)
	simdata1BIN <- summaryBy(PRED~TADBIN, data=simdata1, FUN=median, na.rm=T)
	simdata1BIN

#Merge median PREDs into simulated dataset matching for their TADBIN
	simdata1 <- merge(simdata1,simdata1BIN, by=c("TADBIN"),all=T)
	simdata1 <- rename(simdata1, c("PRED.median" = "PREDMED"))

	simdata1 <- simdata1[with(simdata1, order(simdata1$SIM, simdata1$ID, simdata1$TIME, simdata1$TADBIN)), ]

	orgdata1 <- orgdata1[with(orgdata1, order(orgdata1$ID, orgdata1$TIME, orgdata1$TADBIN)), ]

	orgdata1BIN <- summaryBy(DV~TADBIN, data=orgdata1, FUN=median, na.rm=T)

#Subset for one simulation of the same length of the original dataset
	simdata1ONE <- subset(simdata1, simdata1$SIM == 1)

#Add median PRED for each TADBIN to the orignal dataset
	orgdata1$PREDMED <- simdata1ONE$PREDMED
	orgdata1$PRED <- simdata1ONE$PRED

#-------------------------------------------------------------------------------
#PRED Correction
#Calculate the prediction corrected observed and simulated DVs
	orgdata1$pcY <- (orgdata1$DV)*(orgdata1$PREDMED)/(orgdata1$PRED)
	simdata1$pcY <- (simdata1$DV)*(simdata1$PREDMED)/(simdata1$PRED)

	#Uppsala Style
	#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
	#Plot the confidence interval for the simulated data's percentiles for each bin
	#(for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
	# from all of the simulated datasets compute the 95% CI of these percentiles).

		#Calculate 5, 50 and 95 percentiles for each simulated study (S)
		simdata1.bystudy.median <- ddply(simdata1, .(SIM,TADBIN), function(df) median(df$pcY))
		simdata1.bystudy.median <- rename(simdata1.bystudy.median, c("V1"="medianS"))

		simdata1.bystudy.loCI <- ddply(simdata1, .(SIM,TADBIN), function(df) CI90lo(df$pcY))
		simdata1.bystudy.loCI <- rename(simdata1.bystudy.loCI, c("5%"="loCI90S"))

		simdata1.bystudy.hiCI <- ddply(simdata1, .(SIM,TADBIN), function(df) CI90hi(df$pcY))
		simdata1.bystudy.hiCI <- rename(simdata1.bystudy.hiCI, c("95%"="hiCI90S"))

		simdata1.bystudy <- data.frame(simdata1.bystudy.median, "loCI90S"=simdata1.bystudy.loCI$loCI90S, "hiCI90S"=simdata1.bystudy.hiCI$hiCI90S)

 		titletext <- "PCVPC - Uppsala Style\n"

		plotobj1 <- NULL
		#titletext <- "VPC - Uppsala Style\n"
		plotobj1 <- ggplot(data=orgdata1)
		plotobj1 <- plotobj1 + ggtitle(titletext)

		#Median simulated with confidence band
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata1.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="lightpink")
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata1.bystudy, fun.y=median, geom="line", colour="black", size=1)

		#Lower 90% CI simulated with confidence band
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata1.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata1.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		#Upper 90% CI simulated with confidence band
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata1.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata1.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		plotobj1 <- plotobj1 + geom_point(aes(x=TADBIN, y=pcY), colour="blue", shape = 1)
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=median, geom="line", colour="red", size=1)
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
		plotobj1 <- plotobj1 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
		plotobj1 <- plotobj1 + scale_y_log10("Prediction Corrected (%) \n", lim = c(0.00001, 1000))
		plotobj1 <- plotobj1 + scale_x_continuous("\nTime")
		plotobj1 <- plotobj1 + opts(strip.background = theme_rect(fill = "grey95", colour = "grey50"))
		plotobj1

# ------------------------------------------------------------------------------
# simdata2

#Calculate the median PRED for each TADBIN
	simdata2$PRED <- as.numeric(simdata2$PRED)
	simdata2BIN <- summaryBy(PRED~TADBIN, data=simdata2, FUN=median, na.rm=T)
	simdata2BIN

#Merge median PREDs into simulated dataset matching for their TADBIN
	simdata2 <- merge(simdata2,simdata2BIN, by=c("TADBIN"),all=T)
	simdata2 <- rename(simdata2, c("PRED.median" = "PREDMED"))

	simdata2 <- simdata2[with(simdata2, order(simdata2$SIM, simdata2$ID, simdata2$TIME, simdata2$TADBIN)), ]

	orgdata2 <- orgdata2[with(orgdata2, order(orgdata2$ID, orgdata2$TIME, orgdata2$TADBIN)), ]

	orgdata2BIN <- summaryBy(DV~TADBIN, data=orgdata2, FUN=median, na.rm=T)

#Subset for one simulation of the same length of the original dataset
	simdata2ONE <- subset(simdata2, simdata2$SIM == 1)

#Add median PRED for each TADBIN to the orignal dataset
	orgdata2$PREDMED <- simdata2ONE$PREDMED
	orgdata2$PRED <- simdata2ONE$PRED

#-------------------------------------------------------------------------------
#PRED Correction
#Calculate the prediction corrected observed and simulated DVs
	orgdata2$pcY <- (orgdata2$DV)*(orgdata2$PREDMED)/(orgdata2$PRED)
	simdata2$pcY <- (simdata2$DV)*(simdata2$PREDMED)/(simdata2$PRED)

	#Uppsala Style
	#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
	#Plot the confidence interval for the simulated data's percentiles for each bin
	#(for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
	# from all of the simulated datasets compute the 95% CI of these percentiles).

		#Calculate 5, 50 and 95 percentiles for each simulated study (S)
		simdata2.bystudy.median <- ddply(simdata2, .(SIM,TADBIN), function(df) median(df$pcY))
		simdata2.bystudy.median <- rename(simdata2.bystudy.median, c("V1"="medianS"))

		simdata2.bystudy.loCI <- ddply(simdata2, .(SIM,TADBIN), function(df) CI90lo(df$pcY))
		simdata2.bystudy.loCI <- rename(simdata2.bystudy.loCI, c("5%"="loCI90S"))

		simdata2.bystudy.hiCI <- ddply(simdata2, .(SIM,TADBIN), function(df) CI90hi(df$pcY))
		simdata2.bystudy.hiCI <- rename(simdata2.bystudy.hiCI, c("95%"="hiCI90S"))

		simdata2.bystudy <- data.frame(simdata2.bystudy.median, "loCI90S"=simdata2.bystudy.loCI$loCI90S, "hiCI90S"=simdata2.bystudy.hiCI$hiCI90S)

 		titletext <- "PCVPC - Uppsala Style\n"

		plotobj2 <- NULL
		#titletext <- "VPC - Uppsala Style\n"
		plotobj2 <- ggplot(data=orgdata2)
		plotobj2 <- plotobj2 + ggtitle(titletext)

		#Median simulated with confidence band
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata2.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="lightpink")
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata2.bystudy, fun.y=median, geom="line", colour="black", size=1)

		#Lower 90% CI simulated with confidence band
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata2.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata2.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		#Upper 90% CI simulated with confidence band
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata2.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata2.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		plotobj2 <- plotobj2 + geom_point(aes(x=TADBIN, y=pcY), colour="blue", shape = 1)
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=median, geom="line", colour="red", size=1)
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
		plotobj2 <- plotobj2 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
		plotobj2 <- plotobj2 + scale_y_log10("Prediction Corrected (%) \n", lim = c(0.00001, 1000))
		plotobj2 <- plotobj2 + scale_x_continuous("\nTime")
		plotobj2 <- plotobj2 + opts(strip.background = theme_rect(fill = "grey95", colour = "grey50"))
		plotobj2

# -----------------------------------------------------------------------------
# simdata3

#Calculate the median PRED for each TADBIN
	simdata3$PRED <- as.numeric(simdata3$PRED)
	simdata3BIN <- summaryBy(PRED~TADBIN, data=simdata3, FUN=median, na.rm=T)
	simdata3BIN

#Merge median PREDs into simulated dataset matching for their TADBIN
	simdata3 <- merge(simdata3,simdata3BIN, by=c("TADBIN"),all=T)
	simdata3 <- rename(simdata3, c("PRED.median" = "PREDMED"))

	simdata3 <- simdata3[with(simdata3, order(simdata3$SIM, simdata3$ID, simdata3$TIME, simdata3$TADBIN)), ]

	orgdata3 <- orgdata3[with(orgdata3, order(orgdata3$ID, orgdata3$TIME, orgdata3$TADBIN)), ]

	orgdata3BIN <- summaryBy(DV~TADBIN, data=orgdata3, FUN=median, na.rm=T)

#Subset for one simulation of the same length of the original dataset
	simdata3ONE <- subset(simdata3, simdata3$SIM == 1)

#Add median PRED for each TADBIN to the orignal dataset
	orgdata3$PREDMED <- simdata3ONE$PREDMED
	orgdata3$PRED <- simdata3ONE$PRED

#-------------------------------------------------------------------------------
#PRED Correction
#Calculate the prediction corrected observed and simulated DVs
	orgdata3$pcY <- (orgdata3$DV)*(orgdata3$PREDMED)/(orgdata3$PRED)
	simdata3$pcY <- (simdata3$DV)*(simdata3$PREDMED)/(simdata3$PRED)

	#Uppsala Style
	#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
	#Plot the confidence interval for the simulated data's percentiles for each bin
	#(for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
	# from all of the simulated datasets compute the 95% CI of these percentiles).

		#Calculate 5, 50 and 95 percentiles for each simulated study (S)
		simdata3.bystudy.median <- ddply(simdata3, .(SIM,TADBIN), function(df) median(df$pcY))
		simdata3.bystudy.median <- rename(simdata3.bystudy.median, c("V1"="medianS"))

		simdata3.bystudy.loCI <- ddply(simdata3, .(SIM,TADBIN), function(df) CI90lo(df$pcY))
		simdata3.bystudy.loCI <- rename(simdata3.bystudy.loCI, c("5%"="loCI90S"))

		simdata3.bystudy.hiCI <- ddply(simdata3, .(SIM,TADBIN), function(df) CI90hi(df$pcY))
		simdata3.bystudy.hiCI <- rename(simdata3.bystudy.hiCI, c("95%"="hiCI90S"))

		simdata3.bystudy <- data.frame(simdata3.bystudy.median, "loCI90S"=simdata3.bystudy.loCI$loCI90S, "hiCI90S"=simdata3.bystudy.hiCI$hiCI90S)

 		titletext <- "PCVPC - Uppsala Style\n"

		plotobj3 <- NULL
		#titletext <- "VPC - Uppsala Style\n"
		plotobj3 <- ggplot(data=orgdata3)
		plotobj3 <- plotobj3 + ggtitle(titletext)

		#Median simulated with confidence band
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata3.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="lightpink")
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=medianS), data=simdata3.bystudy, fun.y=median, geom="line", colour="black", size=1)

		#Lower 90% CI simulated with confidence band
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata3.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=loCI90S), data=simdata3.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		#Upper 90% CI simulated with confidence band
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata3.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=hiCI90S), data=simdata3.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

		plotobj3 <- plotobj3 + geom_point(aes(x=TADBIN, y=pcY), colour="blue", shape = 1)
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=median, geom="line", colour="red", size=1)
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
		plotobj3 <- plotobj3 + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
		plotobj3 <- plotobj3 + scale_y_log10("Prediction Corrected (%) \n", lim = c(0.00001, 1000))
		plotobj3 <- plotobj3 + scale_x_continuous("\nTime")
		plotobj3 <- plotobj3 + opts(strip.background = theme_rect(fill = "grey95", colour = "grey50"))
		plotobj3
