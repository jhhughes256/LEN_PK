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

# Extract AUC data
	aucdata <- ddply(fitdata[fitdata$OCC == 12,], .(STUDY, DOSELVL, ID), function(x) {
		if (unique(x$STUDY) %in% c(5115, 6003)) {
			out <- ddply(x[x$TIME == 24, ], .(TIME), function(y) {
		    head(y, n = 1)
		  })
			out$AUC <- out$AUC1
		} else if (unique(x$STUDY) == 8056) {
			out <- ddply(x[x$TIME %in% c(48,72), ], .(TIME), function(y) {
				head(y, n = 2)
			})
			half.nobs <- length(out$AUC2)/2
			out$AUC <- c(out$AUC2[1:half.nobs*2 - 1], out$AUC3[1:half.nobs*2])
		} else {
			out <- ddply(x[x$TIME %in% c(24, 120), ], .(TIME), function(y) {
				head(y, n = 2)
			})
			half.nobs <- length(out$AUC1)/2
			out$AUC <- c(out$AUC1[1:half.nobs*2 - 1], out$AUC5[1:half.nobs*2])
		}
		out
	})

	auc.each <- ddply(aucdata, .(ID), function(x) {
		dose <- if (unique(x$DOSELVL) %in% 1:3) {
			2.5
		} else if (unique(x$DOSELVL) == 4) {
			5
		} else if (unique(x$DOSELVL) == 5) {
			7.5
		} else if (unique(x$DOSELVL) == 6) {
			15
		} else if (unique(x$DOSELVL) == 7) {
			20
		} else if (unique(x$DOSELVL) == 8) {
			25
		} else if (unique(x$DOSELVL) == 9) {
			30
		} else if (unique(x$DOSELVL) == 10) {
			35
		} else if (unique(x$DOSELVL) == 11) {
			50
		} else if (unique(x$DOSELVL) == 12) {
			75
		}
		data.frame(
			ID = unique(x$ID),
			DOSELVL = unique(x$DOSELVL),
			AUC = mean(x$AUC),
			DOSE = dose,
			NORMAUC = mean(x$AUC)/dose
		)
	})

	auc.table <- ddply(auc.each, .(DOSELVL), function(x) {
		data.frame(
			NID = length(x$ID),
			AUC = mean(x$NORMAUC)
		)
	})
