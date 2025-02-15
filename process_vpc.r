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
	master.dir<-"E:/Hughes/Data/PK/FLAG/COV15"
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

#-------------------------------------------------------------------------------
#Visual predictive checks
#-------------------------------------------------------------------------------
#Process the simulated *.fit file.
#Run name - Change this to the RUN you want to process
	runname <- "RUN016_CL_CRCL2_FFM_VPC"

#Process the fit file - Comment this out if you have already generated the csv; this will save time!
  processSIMdata(paste(runname,".ctl",sep=""))    # from the FUNCTION UTILITY

#Read the simulated data
  SIM.data <- read.csv(paste(runname,".nm7/",runname,".fit.csv",sep=""), stringsAsFactors=F, na.strings=".")
  #SIM.data <- rename(SIM.data, c("TAFDE"="TIME"))  # rename time
 length (SIM.data$DV)

 #Subset Simulated data
 SIM.data<- subset(SIM.data, MDV==0)

#-------------------------------------------------------------------------------
#Read the original data
	ORG.data <- read.csv("nmprep_flagged.csv", stringsAsFactors=F, na.strings=".")
  ORG.data <- rename(ORG.data, c("X.ID"="ID"))
  #ORG.data <- rename(ORG.data, c("TAFDE"="TIME"))  # rename time
  ORG.data <- subset (ORG.data, MDV==0 & EVID<=1) #removes data not used in the analysis commented out in the original CTL
	if(str_detect(master.dir, "OSU_SUB")) ORG.data <- subset (ORG.data, STUDY<=6003) # only include if OSU_SUB

 #------------------------------------------------------------------------------
#Assign some factors - ORG.data

	#Bin time
	ORG.data$TIMEBIN <- cut2(ORG.data$TIME, cuts=c(1.02,2.02,4.02,8.02,49,97.5), levels.mean=T)
	ORG.data$TIMEBIN <- as.numeric(paste(ORG.data$TIMEBIN))

	#Bin time after first dose
	ORG.data$TADBIN <- cut2(ORG.data$TAD, cuts=c(1.02,2.02,4.02,8.02,49,97.5), levels.mean=T)
	ORG.data$TADBIN <- as.numeric(paste(ORG.data$TADBIN))

  # Assign some factors
   ORG.data$IDf <- as.factor(ORG.data$ID)
   ORG.data$RACEf <- factor(ORG.data$RACE, labels=c("White", "Other"))
   ORG.data$GENDf <- factor(ORG.data$GEND, labels=c("F", "M"))
   ORG.data$DXCATf <- factor(ORG.data$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
   ORG.data$STUDYf <- factor(ORG.data$STUDY)
	 ORG.data$DOSEf <- factor(ifelse(ORG.data$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	 ORG.data$CRCLf <- factor(ifelse(ORG.data$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	# Assign LLOQ
	 ORG.data$LOQ <- ifelse(ORG.data$STUDY == 6003, 0.005, 0.00025926)

#-------------------------------------------------------------------------------
#Assign some factors - SIM.data

	#Bin time
	SIM.data$TIMEBIN <- cut2(SIM.data$TIME, cuts=c(1.02,2.02,4.02,8.02,49,97.5), levels.mean=T)
	SIM.data$TIMEBIN <- as.numeric(paste(SIM.data$TIMEBIN))

	#Bin time after first dose
	SIM.data$TADBIN <- cut2(SIM.data$TAD, cuts=c(1.02,2.02,4.02,8.02,49,97.5), levels.mean=T)
	SIM.data$TADBIN <- as.numeric(paste(SIM.data$TADBIN))

  # Assign some factors
   SIM.data$IDf <- as.factor(SIM.data$ID)
   SIM.data$RACEf <- factor(SIM.data$RACE, labels=c("White", "Other", "Other"))
   SIM.data$GENDf <- factor(SIM.data$GEND, labels=c("F", "M"))
   SIM.data$DXCATf <- factor(SIM.data$DXCAT, labels=c("CLL","AML","ALL","MM"))
   SIM.data$STUDYf <- factor(SIM.data$STUDY)
	 SIM.data$DOSEf <- factor(ifelse(SIM.data$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))
	 SIM.data$CRCLf <- factor(ifelse(SIM.data$CRCL2<=60,1,2), labels=c("CrCl <60mL/min","CrCl >60mL/min"))

	 # Assign LLOQ
	 SIM.data$LOQ <- ifelse(SIM.data$STUDY == 6003, 0.005, 0.00025926)
	 SIM.data$BLQ <- ifelse(SIM.data$DV < SIM.data$LOQ, 1, 0)

#------------------------------------------------------------------------------
	 #Assign some factors - LOQ.data
		LOQ.data <- read.csv("extval_10156_IBW.csv", stringsAsFactors=F, na.strings=".")
		names(LOQ.data)[1] <- "ID"

	  #Bin time
		LOQ.data <- LOQ.data[LOQ.data$TIME != 0, ]
		LOQ.data$DV[is.na(LOQ.data$DV)] <- 0
	  LOQ.data$TIMEBIN <- cut2(LOQ.data$TIME, cuts=c(1.02,2.02,4.02,8.02,49,97.5), levels.mean=T)
	  LOQ.data$TIMEBIN <- as.numeric(paste(LOQ.data$TIMEBIN))

	  #Bin time after first dose
	  LOQ.data$TADBIN <- cut2(LOQ.data$TAD, cuts=c(0.52,1.02,2.02,4.02,8.02,22.02), levels.mean=T)
	  LOQ.data$TADBIN <- as.numeric(paste(LOQ.data$TADBIN))

	  # Assign some factors
	 	LOQ.data$IDf <- as.factor(LOQ.data$ID)
	 	LOQ.data$RACEf <- factor(LOQ.data$RACE2, labels=c("White", "Other"))
	 	LOQ.data$GENDf <- factor(LOQ.data$GEND, labels=c("F", "M"))
	 	LOQ.data$DXCATf <- factor(LOQ.data$DXCATNUM, labels=c("CLL","AML","ALL","MM"))
	 	LOQ.data$STUDYf <- factor(LOQ.data$STUDY)
	 	LOQ.data$DOSEf <- factor(ifelse(LOQ.data$DOSELVL<=5,1,2), labels=c("Dose <10mg","Dose >10mg"))

		LOQ.data$LOQ <- ifelse(LOQ.data$STUDY == 6003, 0.005, 0.00025926)
		LOQ.data$BLQ <- ifelse(LOQ.data$DV < LOQ.data$LOQ, 1, 0)

#-------------------------------------------------------------------------------
#Confidence intervals - from function utility

	CI90lo <- function(x) quantile(x, probs=0.05)
	CI90hi <- function(x) quantile(x, probs=0.95)

	CI95lo <- function(x) quantile(x, probs=0.025)
	CI95hi <- function(x) quantile(x, probs=0.975)

	#-----------------------------------------------------------------------------

	#Change working directory, so that plots are stored within the VPC folder
	setwd(paste(master.dir,"/",runname,".nm7",sep=""))


#-------------------------------------------------------------------------------
#New style VPC #3 - Uppsala Style
#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
#Plot the confidence interval for the simulated data's percentiles for each bin
# (for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
# from all of the simulated datasets compute the 95% CI of these percentiles).


	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) median(df$DV))
	#SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN, DVIDf), function(df) median(df$DV)) ; DVIDf for metabolite model
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
	#SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN, DVIDf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
	#SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN, DVIDf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

  #Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
	                             atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
	                                  "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
  #plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)


	plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)\n")
  plotobj <- plotobj + scale_x_log10("\nTime after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj

	#normal scale- non-facetted
  ggsave("Uppsala_VPC.png", width=20, height=16, units=c("cm"))

	# log Y-scale-non- facetted
	plotobj2 <- plotobj + scale_y_log10("Concentration (ug/L)\n")
	plotobj2
	ggsave("Uppsala_VPC_log.png", width=20, height=16, units=c("cm"))

	#facet on study---------------------------------------------------------------

	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN,STUDYf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN,STUDYf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN,STUDYf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

  #Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
	                             atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
	                                  "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
  #plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
  plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj <- plotobj + facet_wrap(~STUDYf)
	plotobj
	ggsave("Uppsala_VPC_bystudy.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_bystudylog.png", width=20, height=16, units=c("cm"))

	#facet on dose----------------------------------------------------------------
	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN,DOSEf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN,DOSEf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN,DOSEf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

  #Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
	                             atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
	                                  "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
	  #plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
  plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj <- plotobj + facet_wrap(~DOSEf)
	plotobj
	ggsave("Uppsala_VPC_bydose.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_bydoselog.png", width=20, height=16, units=c("cm"))

	#facet on crcl----------------------------------------------------------------
	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN,CRCLf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN,CRCLf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN,CRCLf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

  #Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
	                             atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
	                                  "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
	  #plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
  plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj <- plotobj + facet_wrap(~CRCLf)
	plotobj
	ggsave("Uppsala_VPC_bycrcl.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_bycrcllog.png", width=20, height=16, units=c("cm"))
#-------------------------------------------------------------------------------
#TADbin VPC

	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN), function(df) median(df$DV))
	#SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN, DVIDf), function(df) median(df$DV)) ; DVIDf for metabolite model
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN), function(df) CI90lo(df$DV))
	#SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN, DVIDf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN), function(df) CI90hi(df$DV))
	#SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN, DVIDf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

  #Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
	                             atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
	                                  "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
  #plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)


	plotobj <- plotobj + geom_point(aes(x=TADBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)\n")
  plotobj <- plotobj + scale_x_log10("\nTime after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj

	#normal scale- non-facetted
  ggsave("Uppsala_VPC_TAD.png", width=20, height=16, units=c("cm"))

	# log Y-scale-non- facetted
	plotobj2 <- NULL
	plotobj2 <- plotobj + scale_y_log10("Concentration (ug/L)\n")
	plotobj2
	ggsave("Uppsala_VPC_TADlog.png", width=20, height=16, units=c("cm"))
	plotobj_VPC <- NULL
	plotobj_VPC <- plotobj2 + scale_x_log10("\nTime after dose (hours)", lim=c(0.3,24))

	plotobj3 <- NULL
	plotobj3 <- plotobj2 + scale_x_continuous("\nTime after dose (hours)")
	plotobj3
	ggsave("Uppsala_VPC_TADlog_nologx.png", width=20, height=16, units=c("cm"))

	#facet on study---------------------------------------------------------------

	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN,STUDYf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN,STUDYf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN,STUDYf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

	#Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
															 atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
																		"Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

	#titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TADBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
	# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	#plotobj <- plotobj + geom_hline(aes(yintercept = LOQ), data = ORG.data, colour = "darkgreen", linetype = "dashed")

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
	plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
	plotobj <- plotobj + facet_wrap(~STUDYf)
	plotobj
	ggsave("Uppsala_VPC_TADbystudy.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_TADbystudylog.png", width=20, height=16, units=c("cm"))

	#facet on dose----------------------------------------------------------------
	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN,DOSEf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN,DOSEf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN,DOSEf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

	#Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
															 atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
																		"Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

	#titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
		#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TADBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
	# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
	plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
	plotobj <- plotobj + facet_wrap(~DOSEf)
	plotobj
	ggsave("Uppsala_VPC_TADbydose.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_TADbydoselog.png", width=20, height=16, units=c("cm"))

	#facet on crcl----------------------------------------------------------------
	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN,CRCLf), function(df) median(df$DV))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN,CRCLf), function(df) CI90lo(df$DV))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN,CRCLf), function(df) CI90hi(df$DV))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

	plotobj <- NULL
	plotobj2 <- NULL

	#Titletext with 3 lines
	titletext <- expression(atop(VPC~plot~Uppsala~Style,
															 atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
																		"Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

	#titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
		#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
	#plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TADBIN, y=DV), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
	# For simulated data: add median, CI90lo, CI90hi
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

	plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
	plotobj <- plotobj + ggtitle(titletext)
	plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)")
	plotobj <- plotobj + scale_x_log10("Time after dose (hours)")
	#plotobj <- plotobj + facet_wrap(~DVIDf)
	plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
	plotobj <- plotobj + facet_wrap(~CRCLf)
	plotobj
	ggsave("Uppsala_VPC_TADbycrcl.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("log (Concentration (ug/L) )")
	plotobj2
	ggsave("Uppsala_VPC_TADbycrcllog.png", width=20, height=16, units=c("cm"))

	#-------------------------------------------------------------------------------
	#TIMEBIN LOQ VPC
	#Set up data for LOQ plots
		LOQ.BLQ.summary <- cbind(
			summaryBy(BLQ~TIMEBIN, data=LOQ.data, FUN=sum),
			summaryBy(DV~TIMEBIN, data=LOQ.data, FUN=length)[-1])
		LOQ.BLQ.summary$LOQPerBLQ <- 100*(LOQ.BLQ.summary$BLQ.sum)/(LOQ.BLQ.summary$DV.length)

	#Simulations calculate mean % BLQ and CI's
		SIM.BLQ.summary <- cbind(
			summaryBy(BLQ~TIMEBIN+SIM, data=SIM.data, FUN=sum),
			summaryBy(DV~TIMEBIN+SIM, data=SIM.data, FUN=length)[-(1:2)])
		SIM.BLQ.summary$PerBLQ <- 100*(SIM.BLQ.summary$BLQ.sum)/(SIM.BLQ.summary$DV.length)

	#Calculated for a table
		SIM.PerBLQ.median <- summaryBy(PerBLQ~TIMEBIN, data=SIM.BLQ.summary, FUN=median)
		SIM.PerBLQ.ci90lo <- summaryBy(PerBLQ~TIMEBIN, data=SIM.BLQ.summary, FUN=CI90lo)
		SIM.PerBLQ.ci90hi <- summaryBy(PerBLQ~TIMEBIN, data=SIM.BLQ.summary, FUN=CI90hi)
		SIM.PerBLQ.summary <- cbind(SIM.PerBLQ.median, SIM.PerBLQ.ci90lo[-1], SIM.PerBLQ.ci90hi[-1])

		ALL.BLQ.summary <- cbind(LOQ.BLQ.summary, SIM.PerBLQ.summary[-1])

#-------------------------------------------------------------------------------
	#Median simulated BLQ with confidence intervals
	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
		SIM.BLQ.summary.bystudy <- ddply(SIM.BLQ.summary, .(SIM,TIMEBIN), function(x) median(x$PerBLQ))
		names(SIM.BLQ.summary.bystudy)[3] <- "medianS"

		plotobj_BLQ <- NULL
		plotobj_BLQ <- ggplot(data=LOQ.BLQ.summary)

	#Median simulated with confidence band
		plotobj_BLQ <- plotobj_BLQ + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.BLQ.summary.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
		plotobj_BLQ <- plotobj_BLQ + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.BLQ.summary.bystudy, geom="line", size = 2, fun.y=median)
		plotobj_BLQ

	#ORG DATA
		plotobj_BLQ <- plotobj_BLQ + geom_point(aes(x=TIMEBIN, y=LOQPerBLQ), shape = 1, colour = "blue", size=2)
		plotobj_BLQ

	#plotobj_BLQ <- plotobj_BLQ + ggtitle(titletext)
    plotobj_BLQ <- plotobj_BLQ + scale_y_continuous("% <LLOQ")
    plotobj_BLQ <- plotobj_BLQ + scale_x_continuous("Time (h)", lim=c(0,160))
	#plotobj_BLQ <- plotobj_BLQ + facet_wrap(~DVIDf)
		plotobj_BLQ <- plotobj_BLQ + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
    plotobj_BLQ

#-------------------------------------------------------------------------------
	#TADBIN LOQ VPC
	#Set up data for LOQ plots
		LOQ.BLQ.summary <- cbind(
			summaryBy(BLQ~TADBIN, data=LOQ.data, FUN=sum),
			summaryBy(DV~TADBIN, data=LOQ.data, FUN=length)[-1])
		LOQ.BLQ.summary$LOQPerBLQ <- 100*(LOQ.BLQ.summary$BLQ.sum)/(LOQ.BLQ.summary$DV.length)

	#Simulations calculate mean % BLQ and CI's
		SIM.BLQ.summary <- cbind(
			summaryBy(BLQ~TADBIN+SIM, data=SIM.data, FUN=sum),
			summaryBy(DV~TADBIN+SIM, data=SIM.data, FUN=length)[-(1:2)])
		SIM.BLQ.summary$PerBLQ <- 100*(SIM.BLQ.summary$BLQ.sum)/(SIM.BLQ.summary$DV.length)

	#Calculated for a table
		SIM.PerBLQ.median <- summaryBy(PerBLQ~TADBIN, data=SIM.BLQ.summary, FUN=median)
		SIM.PerBLQ.ci90lo <- summaryBy(PerBLQ~TADBIN, data=SIM.BLQ.summary, FUN=CI90lo)
		SIM.PerBLQ.ci90hi <- summaryBy(PerBLQ~TADBIN, data=SIM.BLQ.summary, FUN=CI90hi)
		SIM.PerBLQ.summary <- cbind(SIM.PerBLQ.median, SIM.PerBLQ.ci90lo[-1], SIM.PerBLQ.ci90hi[-1])

		ALL.BLQ.summary <- cbind(LOQ.BLQ.summary, SIM.PerBLQ.summary[-1])

#-------------------------------------------------------------------------------
	#TADBIN
		SIM.BLQ.summary.bystudy <- ddply(SIM.BLQ.summary, .(SIM,TADBIN), function(x) median(x$PerBLQ))
		names(SIM.BLQ.summary.bystudy)[3] <- "medianS"

		plotobj_BLQ <- NULL
		plotobj_BLQ <- ggplot(data=LOQ.BLQ.summary)

	#Median simulated with confidence band
		plotobj_BLQ <- plotobj_BLQ + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.BLQ.summary.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
		plotobj_BLQ <- plotobj_BLQ + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.BLQ.summary.bystudy, geom="line", size = 2, fun.y=median)
		plotobj_BLQ

	#ORG DATA
		plotobj_BLQ <- plotobj_BLQ + geom_point(aes(x=TADBIN, y=LOQPerBLQ), shape = 1, colour = "blue", size=2)
		plotobj_BLQ

	#plotobj_BLQ <- plotobj_BLQ + ggtitle(titletext)
    plotobj_BLQ <- plotobj_BLQ + scale_y_continuous("% <LLOQ\n\n")
    plotobj_BLQ <- plotobj_BLQ + scale_x_log10("\nTime (h)", lim=c(0.3,24))
	#plotobj_BLQ <- plotobj_BLQ + facet_wrap(~DVIDf)
		plotobj_BLQ <- plotobj_BLQ + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
    plotobj_BLQ

		ggsave("Percent_BLOQ_Plot.png", width=10, height=5, units=c("cm"))

#-------------------------------------------------------------------------------
	#Stack the BLQ plot and VCP plots
	filename <- "VPC_BLQ_Stack_study.png"
	png(filename, width = 720, height = 920)

	# 2 ggplot2 graphs in a grid layout
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(3,1)))

  #Plot 1
  plotobj_stack1 <- NULL
  plotobj_stack2 <- NULL
  plotobj_stack1 <- plotobj_BLQ
  print(plotobj_stack1, vp=vplayout(3,1))

  #Plot 2
  plotobj_stack2 <- plotobj_VPC
  print(plotobj_stack2, vp=vplayout(1:2,1))

	dev.off()

#-------------------------------------------------------------------------------
###GENERATE A pcVPC
###Bergstrand et al 2011 - Prediction-Corrected Visual Predictive Checks for Diagnosing Nonlinear Mixed-Effects Models

#Calculate the median PRED for each TIMEBIN
	SIM.data$PRED <- as.numeric(SIM.data$PRED)
	SIM.dataBIN <- summaryBy(PRED~TIMEBIN, data=SIM.data, FUN=median, na.rm=T)
	SIM.dataBIN

#Merge median PREDs into simulated dataset matching for their TIMEBIN
	SIM.data <- merge(SIM.data,SIM.dataBIN, by=c("TIMEBIN"),all=T)
	SIM.data <- rename(SIM.data, c("PRED.median" = "PREDMED"))

	SIM.data <- SIM.data[with(SIM.data, order(SIM.data$SIM, SIM.data$ID, SIM.data$TIME, SIM.data$TIMEBIN)), ]

	ORG.data <- ORG.data[with(ORG.data, order(ORG.data$ID, ORG.data$TIME, ORG.data$TIMEBIN)), ]

	ORG.dataBIN <- summaryBy(DV~TIMEBIN, data=ORG.data, FUN=median, na.rm=T)

#Subset for one simulation of the same length of the original dataset
	SIM.dataONE <- subset(SIM.data, SIM.data$SIM == 1)

#Add median PRED for each TIMEBIN to the orignal dataset
	ORG.data$PREDMED <- SIM.dataONE$PREDMED
	ORG.data$PRED <- SIM.dataONE$PRED

#-------------------------------------------------------------------------------
#PRED Correction
#Calculate the prediction corrected observed and simulated DVs
	ORG.data$pcY <- (ORG.data$DV)*(ORG.data$PREDMED)/(ORG.data$PRED)
	SIM.data$pcY <- (SIM.data$DV)*(SIM.data$PREDMED)/(SIM.data$PRED)

#Uppsala Style
#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
#Plot the confidence interval for the simulated data's percentiles for each bin
#(for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
# from all of the simulated datasets compute the 95% CI of these percentiles).

	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) median(df$pcY))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90lo(df$pcY))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90hi(df$pcY))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

		titletext <- "PCVPC - Uppsala Style\n"

	plotobj <- NULL
	#titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)
	plotobj <- plotobj + ggtitle(titletext)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="lightpink")
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=pcY), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=pcY), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=pcY), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=pcY), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + scale_y_log10("Prediction Corrected (%) \n")
	plotobj <- plotobj + scale_x_log10("\nTime")
	plotobj <- plotobj + opts(strip.background = theme_rect(fill = "grey95", colour = "grey50"))
	plotobj

	ggsave("Uppsala_PCVPC_log.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("Prediction Corrected (%) \n", limits = c(0.0001,10))
	plotobj2

	ggsave("Uppsala_PCVPC_adjlog.png", width=20, height=16, units=c("cm"))

	plotobj3 <- plotobj2 + scale_x_continuous("\nTime")
	plotobj3

	ggsave("Uppsala_PCVPC_adjlog_contx.png", width=20, height=16, units=c("cm"))

#-------------------------------------------------------------------------------
###GENERATE A pcVPC
###Bergstrand et al 2011 - Prediction-Corrected Visual Predictive Checks for Diagnosing Nonlinear Mixed-Effects Models

#Calculate the median PRED for each TADBIN
	SIM.dataBIN <- summaryBy(PRED~TADBIN, data=SIM.data, FUN=median, na.rm=T)
	SIM.dataBIN

	SIM.data <- SIM.data[-c(38:39)]

#Merge median PREDs into simulated dataset matching for their TADBIN
	SIM.data <- merge(SIM.data,SIM.dataBIN, by=c("TADBIN"),all=T)
	SIM.data <- rename(SIM.data, c("PRED.median" = "PREDMED"))

	SIM.data <- SIM.data[with(SIM.data, order(SIM.data$SIM, SIM.data$ID, SIM.data$TIME, SIM.data$TADBIN)), ]

	ORG.data <- ORG.data[with(ORG.data, order(ORG.data$ID, ORG.data$TIME, ORG.data$TADBIN)), ]

	ORG.dataBIN <- summaryBy(DV~TADBIN, data=ORG.data, FUN=median, na.rm=T)

#Subset for one simulation of the same length of the original dataset
	SIM.dataONE <- subset(SIM.data, SIM.data$SIM == 1)

#Add median PRED for each TADBIN to the orignal dataset
	ORG.data$PREDMED <- SIM.dataONE$PREDMED
	ORG.data$PRED <- SIM.dataONE$PRED

#-------------------------------------------------------------------------------
#PRED Correction
#Calculate the prediction corrected observed and simulated DVs
	ORG.data$pcY <- (ORG.data$DV)*(ORG.data$PREDMED)/(ORG.data$PRED)
	SIM.data$pcY <- (SIM.data$DV)*(SIM.data$PREDMED)/(SIM.data$PRED)

#Uppsala Style
#Xpose method - http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
#Plot the confidence interval for the simulated data's percentiles for each bin
#(for each simulated data set compute the percentiles for each bin, then, from all of the percentiles
# from all of the simulated datasets compute the 95% CI of these percentiles).

	#Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TADBIN), function(df) median(df$pcY))
	SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

	SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TADBIN), function(df) CI90lo(df$pcY))
	SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

	SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TADBIN), function(df) CI90hi(df$pcY))
	SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

	SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

		titletext <- "PCVPC - Uppsala Style\n"

	plotobj <- NULL
	#titletext <- "VPC - Uppsala Style\n"
	plotobj <- ggplot(data=ORG.data)
	plotobj <- plotobj + ggtitle(titletext)

	#Median simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="lightpink")
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

	#Lower 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	#Upper 90% CI simulated with confidence band
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", fill="skyblue1")
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

	plotobj <- plotobj + geom_point(aes(x=TADBIN, y=pcY), colour="blue", shape = 1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=pcY), fun.y=median, geom="line", colour="red", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + stat_summary(aes(x=TADBIN, y=pcY), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
	plotobj <- plotobj + scale_y_log10("Prediction Corrected (%) \n")
	plotobj <- plotobj + scale_x_log10("\nTime")
	plotobj <- plotobj + opts(strip.background = theme_rect(fill = "grey95", colour = "grey50"))
	plotobj

	ggsave("Uppsala_PCVPC_TADlog.png", width=20, height=16, units=c("cm"))

	plotobj2 <- plotobj + scale_y_log10("Prediction Corrected (%) \n", limits = c(0.0001,10))
	plotobj2

	ggsave("Uppsala_PCVPC_TADadjlog.png", width=20, height=16, units=c("cm"))

	plotobj3 <- plotobj2 + scale_x_continuous("\nTime")
	plotobj3

	ggsave("Uppsala_PCVPC_TADadjlog_contx.png", width=20, height=16, units=c("cm"))

#-------------------------------------------------------------------------------
### NPDE - Normalised Prediction Distribution Errors
### Who knows how this actually works

	ORG.data.NPDE <- ORG.data[c("ID", "TIME", "DV", "MDV")]
	ORG.data.NPDE <- ORG.data.NPDE[ORG.data.NPDE$MDV==0, ]
	ORG.data.NPDE <- data.frame("id"=ORG.data.NPDE$ID,
	                            "xobs"=ORG.data.NPDE$TIME,
	                            "yobs"=ORG.data.NPDE$DV)
	#write.table(ORG.data.NPDE, "ORG.data.NPDE.txt", sep="\t", row.names=F)

	SIM.data.NPDE <- SIM.data[c("ID", "TIME", "DV", "MDV")]
	SIM.data.NPDE <- SIM.data.NPDE[SIM.data.NPDE$MDV==0, ]
	SIM.data.NPDE <- data.frame("idsim"=SIM.data.NPDE$ID,
	                            "xsim"=SIM.data.NPDE$TIME,
	                            "ysim"=SIM.data.NPDE$DV)
	#write.table(SIM.data.NPDE, "SIM.data.NPDE.txt", sep="\t", row.names=F)

	x <-autonpde(ORG.data.NPDE,SIM.data.NPDE,iid=1,ix=2,iy=3,
	                  namsav="npde_results", units=list(x="hr",y="ug/L"), decorr.method="inverse")

	#Save the plots in a picture with all 4 plots above each other-----------------

	#Calculate Normalised Prediction Distribution Errors

		# head(x["results"]["res"])
					# ypred ycomp    pd     ydobs      npde
		# 1 0.7099454  1.00 0.710 0.5235378 0.5533847
		# 2 0.6375447  1.00 0.750 0.3774446 0.5533847
		# 3 0.6092997  1.00 0.775 0.2431595 0.4676988
		# 4 0.5845265  1.00 0.775 0.2274306 0.4124631
		# 5 0.5224427  1.38 0.890 1.6177910 1.5141019
		# 6 0.4136062  1.13 0.900 0.8453862 0.9944579

		#head(x["results"]["res"])$npde

		plot(x)
		y <- gof.test(x)
		out <- capture.output(y)
		cat(out,file = paste(runname, "gof_test.txt", sep = "_"), sep="\n",append=TRUE)

		#filename <- paste(runname,"NPDE.pdf",sep="_")
		#pdf(filename, width = 12, height = 8, pointsize = 10)
		filename <- paste(runname,"NPDE.png",sep="_")
		png(filename, width = 12, height = 12, units = "in", res = 400, pointsize = 10)

		#4 ggplot2 graphs in a grid layout
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(8,12)))

	#QQ plot of NPDEs
	plotobj1 <- NULL
	plotobj1 <- ggplot(x["results"]["res"])
	plotobj1 <- plotobj1 + ggtitle("(a)\n")
	plotobj1 <- plotobj1 + stat_qq(aes(sample = npde), geom = "point", colour = "blue", shape = 1)
	plotobj1 <- plotobj1 + geom_abline(aes(sample = npde), intercept = 0, slope = 1)
	#plotobj1 <- plotobj1 + geom_hline(aes(yintercept = -1.96), linetype = "dashed", alpha = 0.2)
	#plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 1.96), linetype = "dashed", alpha = 0.2)
	#plotobj1 <- plotobj1 + geom_vline(aes(xintercept = -1.96), linetype = "dashed", alpha = 0.2)
	#plotobj1 <- plotobj1 + geom_vline(aes(xintercept = 1.96), linetype = "dashed", alpha = 0.2)
	plotobj1 <- plotobj1 + scale_x_continuous("\nTheoretical Quantiles", lim = c(-3,3))
	plotobj1 <- plotobj1 + scale_y_continuous("Observed Quantiles (NPDE)\n", lim = c(-3,3))
	print(plotobj1, vp=vplayout(1:4,1:6))

	#Distribution of NPDEs versus normal distribution
	norm.dist <- rnorm(n = length(x["results"]["res"])*1000, mean = 0, sd = 1)
	norm.dist <- as.data.frame(norm.dist)

	plotobj2 <- NULL
	plotobj2 <- ggplot(x["results"]["res"])
	plotobj2 <- plotobj2 + ggtitle("(b)\n")
	plotobj2 <- plotobj2 + geom_density(aes(x = npde, y = ..density..), colour = "blue")
	plotobj2 <- plotobj2 + geom_density(aes(x = norm.dist, y = ..density..),data = norm.dist, colour = "red")
	plotobj2 <- plotobj2 + scale_x_continuous("\nNPDE", lim = c(-3,3))
	plotobj2 <- plotobj2 + scale_y_continuous("Distribution Density\n")
	print(plotobj2, vp=vplayout(1:4,7:12))

	#NPDE versus TIME
	x["results"]["res"]$time <- ORG.data.NPDE$xobs

	plotobj3 <- NULL
	plotobj3 <- ggplot(x["results"]["res"])
	plotobj3 <- plotobj3 + ggtitle("(c)\n")
	plotobj3 <- plotobj3 + geom_point(aes(x = time, y = npde), colour = "blue", shape = 1)
	plotobj3 <- plotobj3 + geom_abline(aes(intercept = 0, slope = 0))
	#plotobj3 <- plotobj3 + geom_abline(aes(intercept = -1.96, slope = 0), linetype = "dashed", alpha = 0.2)
	#plotobj3 <- plotobj3 + geom_abline(aes(intercept = 1.96, slope = 0), linetype = "dashed", alpha = 0.2)
	plotobj3 <- plotobj3 + geom_smooth(aes(x = time, y = npde), method = loess, se = F, colour = "red")
	plotobj3 <- plotobj3 + scale_x_continuous("\nHours Since First Dose")
	plotobj3 <- plotobj3 + scale_y_continuous("NPDE\n")
	print(plotobj3, vp=vplayout(5:8,1:6))

	#NPDE versus Predictions (ypred)
	plotobj4 <- NULL
	plotobj4 <- ggplot(x["results"]["res"])
	plotobj4 <- plotobj4 + ggtitle("(d)\n")
	plotobj4 <- plotobj4 + geom_point(aes(x = ypred, y = npde), colour = "blue", shape = 1)
	plotobj4 <- plotobj4 + geom_abline(aes(intercept = 0, slope = 0))
	#plotobj4 <- plotobj4 + geom_abline(aes(intercept = -1.96, slope = 0), linetype = "dashed", alpha = 0.2)
	#plotobj4 <- plotobj4 + geom_abline(aes(intercept = 1.96, slope = 0), linetype = "dashed", alpha = 0.2)
	plotobj4 <- plotobj4 + geom_smooth(aes(x = ypred, y = npde), method = loess, se = F, colour = "red")
	plotobj4 <- plotobj4 + scale_x_continuous("\nPredicted (Sample Mean)")
	plotobj4 <- plotobj4 + scale_y_continuous("NPDE\n")
	print(plotobj4, vp=vplayout(5:8,7:12))

	dev.off()

	#Original plots
	#png("npdeplots.png", width = 920, height = 920)
	#par(mfrow=c(2,2))
	#plot(x,plot.type="qqplot", new=F)
	#plot(x,plot.type="hist", new=F)
	#plot(x,plot.type="x.scatter", new=F)
	#plot(x,plot.type="pred.scatter", new=F)
	#dev.off()

	#-------------------------------------------------------------------------------
	### LLOQ VPC
	### To determine the density of LOQ values across the board
	LOQ.data$LOQ <- ifelse(LOQ.data$STUDY == 6003, 0.005, 0.00025926)
	LOQ.data$DV[is.na(LOQ.data$DV)] <- 0

	#Bin time
	LOQ.data$TIMEBIN <- cut2(LOQ.data$TIME, g=10, levels.mean=T)
	LOQ.data$TIMEBIN <- as.numeric(paste(LOQ.data$TIMEBIN))

	#Bin time after first dose
	LOQ.data$TADBIN <- cut2(LOQ.data$TAD, g=10, levels.mean=T)
	LOQ.data$TADBIN <- as.numeric(paste(LOQ.data$TADBIN))

	LOQ.data$BLQ <- ifelse(LOQ.data$DV < LOQ.data$LOQ, 1, 0)
	blqtime <- ldply(unique(LOQ.data$TIMEBIN), function(x) {
		sub <- LOQ.data[LOQ.data$TIMEBIN == x, ]
		data.frame(
			"BLQ" = mean(sub$BLQ),
			"TIMEBIN" = x)
	})
	blqtad <- ldply(unique(LOQ.data$TADBIN), function(x) {
		sub <- LOQ.data[LOQ.data$TADBIN == x, ]
		data.frame(
			"BLQ" = mean(sub$BLQ),
			"TADBIN" = x)
	})

	plotobj2 <- NULL
	plotobj2 <- ggplot(blqtad, aes(TADBIN, BLQ))
	plotobj2 <- plotobj2 + ggtitle("LOQ Plot\n")
	plotobj2 <- plotobj2 + geom_bar(stat = "identity", width = 0.5)
	plotobj2
  ggsave("loq_barplot.png", width=20, height=16, units=c("cm"))
