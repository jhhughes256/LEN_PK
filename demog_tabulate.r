###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "demog_tabulate"
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
  study.numbers <- c("05115", "06003", "08056", "10016", "10156")
  file.dir <- paste0(working.dir, "/datacheck_clin_", study.numbers, "_Output")
  file.names <- paste0(file.dir, "/", study.numbers, "_covdata.csv")

  data <- llply(file.names, function(x) {
    read.csv(x, stringsAsFactors = F, na.strings = c("."))
  })

  data <- llply(data, function(x) {
    if (length(x$SECRMGDL) == 0) {
      x$CRCL <- 0
    } else {
      fsex <- ifelse(x$GEND==1,1.23,1.04)
      ibw <- ifelse(x$GEND==1,50+0.9*(x$HT-152),45.5+0.9*(x$HT-152))
      x$CRCL <- (140-x$AGE)*ibw*fsex/(x$SECRMGDL*88.4)
    }
    x
  })

  ldply(data, function(x) {
    c(no = length(x$ID), nsex = length(x$GEND[x$GEND == 0]), msex = 1 - mean(x$GEND),
      minage = min(x$AGE), medage = median(x$AGE), maxage = max(x$AGE),
      minwt = min(x$WT), medwt = median(x$WT), maxwt = max(x$WT),
      mincrcl = min(x$CRCL), medcrcl = median(x$CRCL), maxcrcl = max(x$CRCL),
      mindose = min(x$DOSEMG), meddose = median(x$DOSEMG), maxdose = max(x$DOSEMG)
    )
  })

## -----------------------------------------------------------------------------
  all.studies <- rbind(data[[1]], data[[2]], data[[3]], data[[4]])
  sum.fun.all <- function(x) {
    c(no = length(x$ID), nsex = length(x$GEND[x$GEND == 0]), msex = 1 - mean(x$GEND),
      minage = min(x$AGE), medage = median(x$AGE), maxage = max(x$AGE),
      minwt = min(x$WT), medwt = median(x$WT), maxwt = max(x$WT),
      mincrcl = min(x$CRCL), medcrcl = median(x$CRCL), maxcrcl = max(x$CRCL),
      mindose = min(x$DOSEMG), meddose = median(x$DOSEMG), maxdose = max(x$DOSEMG)
    )
  }
  sum.fun.all(all.studies)
  with(all.studies, table(DXCATNUM))
