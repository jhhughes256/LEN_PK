###process_fit.r
##Goal: To process fit file form nonmem for analysis of model
##Note: Based heavily off of process_fit_10July16.r -> DF code

#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()

# Set the working directory
   master.dir <- "E:/Hughes/Data/PK/REDO"
   scriptname <- "process_fit"
   setwd(master.dir)

 #Load libraries
  library(ggplot2)
  library(doBy)
  library(plyr)
  library(stringr)
  library(R2HTML)
  library(Hmisc)
  library(grid)
  library(reshape2)
  library(GGally)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

#Use custom ggplot 2 theme
	theme_custom <- theme_set(theme_bw(18))
	theme_custom <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"))

#-------------------------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
	theme_bw2 <- theme_set(theme_bw(base_size = 22))
	theme_bw2 <- theme_update(plot.margin = unit(c(1,1,3,1), "lines"),
														axis.title.x=element_text(size = 16, vjust = 0),
														axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
														strip.text.x=element_text(size = 16),
														strip.text.y=element_text(size = 16, angle = 90))

#-------------------------------------------------------------------------------

#Define metadata  #GOLD

#Comment
	metacomment <- NULL
	metacomment <-  "Lenalidomide PK"

#Generating script - this is added using a macro defined in Notepad++
	metafilepath <- NULL
	metafilepath <- "E:/Hughes/lena_scripts/process_fit.r"
	scriptname <- gsub(".r","",basename(metafilepath))
	working.dir <- dirname(metafilepath)

#Last run-time of script
	metadatetime <- NULL
	metadatetime <- Sys.time()

	metadata1 <- paste("#Comment: ",metacomment,sep="")
	metadata2 <- paste("#Generating script: ",metafilepath,sep="")
	metadata3 <- paste("#Script run: ",metadatetime,sep="")
	metadata4 <- paste("#Script run: ",metadatetime,sep="")
#-------------------------------------------------------------------------------

#Read in the fit file
#Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in,".csv", sep="")

	runfolder <- base.path #picks up the folder of the run being analysed

#Read in process_fit.r & nmprep_clin.r output
  fitdata <- read.csv(file.name.out)[-1]
	nmprep <- read.csv("E:/Hughes/Data/PK/REDO/COV24/nmprep_allstudies.csv")

#Remove dose events & missing values
  fitdata <- subset(fitdata, MDV==0)
	fitdata$HT[fitdata$HT==0&fitdata$GEND==1] <- 1.75
	fitdata$HT[fitdata$HT==0&fitdata$GEND==0] <- 1.6

#Set factors and levels (not all categorical covariates have been included!)
	fitdata$IDf <- as.factor(fitdata$ID)

	fitdata$STUDYf <- factor(fitdata$STUDY)
	levels(fitdata$STUDYf) <- paste("Study",levels(fitdata$STUDYf))

	fitdata$OCCf <- factor(nmprep$OCC[nmprep$MDV == 0])
	levels(fitdata$OCCf) <- paste("Occasion",levels(fitdata$OCCf))

	fitdata$GRPf <- factor(fitdata$GRP)
	levels(fitdata$GRPf) <- paste("Group",levels(fitdata$GRPf))

	fitdata$DOSELVLf <- factor(fitdata$DOSELVL)
	levels(fitdata$DOSELVLf) <- paste(
		c(2.5,2.5,2.5,5,15,20,25,30,35,50,75),
		c("mg QD","mg-5mg QD","mg-7.5mg QD",rep("mg QD",8)),
		sep="")
	fitdata$GENDf <- factor(fitdata$GEND)
	levels(fitdata$GENDf) <- c("F","M")

	fitdata$RACEf <- factor(fitdata$RACE)
	levels(fitdata$RACEf) <- c("Unknown","Caucasian","Non-Caucasian")

	fitdata$DXCATf <- factor(fitdata$DXCAT)
	levels(fitdata$DXCATf) <- c("CLL","AML","ALL","MM")

  fitdataone <- ddply(fitdata, .(ID), oneperID)

#-------------------------------------------------------------------------------
# 10016 Genetic comparisons
#Read in pd data from 10016
	datapd <- read.csv("E:/Hughes/Data/RAW_Clinical/datacheck_pd_10016_Output/10016_allpd.csv")[-1]
	names(datapd)[1] <- "ID"
	datapk <- fitdataone[fitdataone$STUDY == 10016, ]
	datapk$ID <- datapk$ID - 121
	dataplot <- merge(datapd, datapk[-(2:11)])
	sub <- dataplot[dataplot$TRT == "UT" & dataplot$Day == 5,]

	titletext <- "Previously Treated Bone Marrow"
	plotobj1 <- NULL
	plotobj1 <- ggplot(sub)
	plotobj1 <- plotobj1 + ggtitle(titletext)
	plotobj1 <- plotobj1 + geom_point(aes(ddCt, CL), na.rm=TRUE)
	plotobj1 <- plotobj1 + scale_x_continuous(name="ddCt")
	plotobj1 <- plotobj1 + scale_y_continuous(name="CL")
	#plotobj1 <- plotobj1 +
	#plotobj1 <- plotobj1 +
	plotobj1 <- plotobj1 + facet_wrap(~Gene)
	plotobj1
