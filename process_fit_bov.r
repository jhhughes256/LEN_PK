###process_fit_bov.r
##Goal: To process fit file for nonmem for analysis of model specfically for BOV

#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/REDO"
  scriptname <- "process_fit_bov"
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
#stops the y axis lable mashing onto the axis values			#NOTE: untested code from DF
  #theme_custom <- theme_update(axis.title.y=element_text(vjust=1, angle = 90))

#--------------------------------------------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
	theme_bw2 <- theme_set(theme_bw(base_size = 22))
	theme_bw2 <- theme_update(plot.margin = unit(c(1,1,3,1), "lines"),
														axis.title.x=element_text(size = 16, vjust = 0),
														axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
														strip.text.x=element_text(size = 16),
														strip.text.y=element_text(size = 16, angle = 90))

#--------------------------------------------------------------------------------------------------

#Define metadata  #GOLD

#Comment
	metacomment <- NULL
	metacomment <-  "Lenalidomide PK"

#Generating script - this is added using a macro defined in Notepad++ - invoked using ctl-alt-M
	metafilepath <- NULL
	metafilepath <- "E:/Hughes/lena_scripts/process_fit_bov.r"
	scriptname <- gsub(".r","",basename(metafilepath))
	working.dir <- dirname(metafilepath)

#Last run-time of script
	metadatetime <- NULL
	metadatetime <- Sys.time()

	metadata1 <- paste("#Comment: ",metacomment,sep="")
	metadata2 <- paste("#Generating script: ",metafilepath,sep="")
	metadata3 <- paste("#Script run: ",metadatetime,sep="")
	metadata4 <- paste("#Script run: ",metadatetime,sep="")
#--------------------------------------------------------------------------------------------------

#Read in the fit file
#Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in,".csv", sep="")

	runfolder <- base.path #picks up the folder of the run being analysed


#Read *.fit file and attach, so column names are available
  fitdata <- read.table(file=file.name.in, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
#Write to file
  write.csv(fitdata, file=file.name.out)

#Remove dose events & missing values
  fitdata <- subset(fitdata, MDV==0)
	fitdata$HT[fitdata$HT==0&fitdata$GEND==1] <- 1.75
	fitdata$HT[fitdata$HT==0&fitdata$GEND==0] <- 1.6

#Set factors and levels (not all categorical covariates have been included!)
	fitdata$IDf <- as.factor(fitdata$ID)

	fitdata$STUDYf <- factor(fitdata$STUDY)
	levels(fitdata$STUDYf) <- paste("Study",levels(fitdata$STUDYf))

	fitdata$OCCf <- factor(fitdata$OCC)
	levels(fitdata$OCCf) <- paste("Occasion",levels(fitdata$OCCf))

	fitdata$GRPf <- factor(fitdata$GRP)
	levels(fitdata$GRPf) <- paste("Group",levels(fitdata$GRPf))

	fitdata$DOSELVLf <- factor(fitdata$DOSELVL)
	levels(fitdata$DOSELVLf) <- paste(c(2.5,2.5,2.5,5,15,20,25,30,35,50,75),
																		c("mg QD","mg-5mg QD","mg-7.5mg QD",rep("mg QD",8)),sep="")
	fitdata$GENDf <- factor(fitdata$GEND)
	levels(fitdata$GENDf) <- c("F","M")

	fitdata$RACEf <- factor(fitdata$RACE)
	levels(fitdata$RACEf) <- c("Unknown","Caucasian","Non-Caucasian")

	fitdata$DXCATf <- factor(fitdata$DXCAT)
	levels(fitdata$DXCATf) <- c("CLL","AML","ALL","MM")

  fitdataone <- ddply(fitdata, .(ID, OCC), oneperID)
#--------------------------------------------------------------------------------------------------
  plotobj <- NULL
  plotobj <- ggplot(data=fitdataone)
  plotobj <- plotobj + geom_boxplot(aes(x = DOSELVLf, y = CL), position=position_dodge(width=0.9))
  plotobj <- plotobj + stat_summary(aes(x = DOSELVLf, y = CL), fun.data = boxplot.give.n, geom = "text", size = 6, colour="red")
  plotobj <- plotobj + scale_x_discrete("Occasions")
  plotobj <- plotobj + scale_y_continuous("ETA4")
  #plotobj  <- plotobj + ggtitle("Final PK model\n")  #legend.position="none",
  #plotobj <- plotobj + ggtitle("Base PK model\n")  #legend.position="none",
  plotobj <- plotobj + ggtitle("ETA PLOT\n")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  plotobj
