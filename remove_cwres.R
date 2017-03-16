###remove_cwres.r
##Goal: To create a file for fitting without high CWRES values

#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()

# Set the working directory
   master.dir <- "E:/Hughes/Data/PK/REDO"
   scriptname <- "remove_cwres"
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
	metacomment <- "Lenalidomide PK"

#Generating script - this is added using a macro defined in Notepad++ - invoked using ctl-alt-M
	metafilepath <- NULL
	metafilepath <- "E:/Hughes/lena_scripts/remove_cwres.r"
	scriptname <- gsub(".r","",basename(metafilepath))
	working.dir <- dirname(metafilepath)

#Last run-time of script
	metadatetime <- NULL
	metadatetime <- Sys.time()

	metadata1 <- paste("#Comment: ",metacomment,sep="")
	metadata2 <- paste("#Generating script: ",metafilepath,sep="")
	metadata3 <- paste("#Script run: ",metadatetime,sep="")
	metadata4 <- paste("#Script run: ",metadatetime,sep="")

# ------------------------------------------------------------------------------
#Read in the fit file
#Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in,".csv", sep="")

  fitdata <- read.csv(file.name.out)
  nmprep <- read.csv("E:/Hughes/Data/PK/REDO/COV24/nmprep_allstudies.csv")

  new.nmprep <- nmprep[-which(fitdata$CWRES > 3.5), ]
  names(new.nmprep)[1] <- "#ID"
  write.csv(new.nmprep, "nmprep_cwresrun058.csv", quote = F, row.names = F)
