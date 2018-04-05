###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "nmprep_earlymdv"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)

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

# Load data
  filename.org <- "E:/Hughes/Data/RAW_Clinical/nmprep_clin_Output/nmprep_flagged_BLQ.csv"
  filename.all <- "E:/Hughes/Data/RAW_Clinical/nmprep_newclin_Output/allnmprep_flagged_BLQ.csv"
  orgdata <- read.csv(filename.org, stringsAsFactors = F)
  alldata <- read.csv(filename.all, stringsAsFactors = F)

# Set MDV before cmax to lloq or half of lloq
  orgdata$OCC <- ceiling(orgdata$TIME/24 + 0.00000001)
  alldata$OCC <- ceiling(alldata$TIME/24 + 0.00000001)
  orgdataLOQ <- ddply(orgdata, .(ID, OCC), function(x) {
    tmax <- x$TAD[which(x$DV == max(x$DV))]
    whichdv <- which(x$AMT == "." & x$TAD < tmax & x$DV == ".")
    if (unique(x$STUDY) == 6003) {
      x$DV[whichdv] <- 5e-3
    } else {
      x$DV[whichdv] <- 2.5926e-4
    }
    x$MDV[whichdv] <- 0
    x
  })
  orgdataHLOQ <- ddply(orgdata, .(ID, OCC), function(x) {
    tmax <- x$TAD[which(x$DV == max(x$DV))]
    whichdv <- which(x$AMT == "." & x$TAD < tmax & x$DV == ".")
    if (unique(x$STUDY) == 6003) {
      x$DV[whichdv] <- 5e-3/2
    } else {
      x$DV[whichdv] <- 2.5926e-4/2
    }
    x$MDV[whichdv] <- 0
    x
  })
  alldataLOQ <- ddply(alldata, .(ID, OCC), function(x) {
    tmax <- x$TAD[which(x$DV == max(x$DV))]
    whichdv <- which(x$AMT == "." & x$TAD < tmax & x$DV == ".")
    if (unique(x$STUDY) == 6003) {
      x$DV[whichdv] <- 5e-3
    } else {
      x$DV[whichdv] <- 2.5926e-4
    }
    x$MDV[whichdv] <- 0
    x
  })
  alldataHLOQ <- ddply(alldata, .(ID, OCC), function(x) {
    tmax <- x$TAD[which(x$DV == max(x$DV))]
    whichdv <- which(x$AMT == "." & x$TAD < tmax & x$DV == ".")
    if (unique(x$STUDY) == 6003) {
      x$DV[whichdv] <- 5e-3/2
    } else {
      x$DV[whichdv] <- 2.5926e-4/2
    }
    x$MDV[whichdv] <- 0
    x
  })

# Save datasets
  filename.org <- paste(output.dir,"nmprep_earlymdv_",sep="/")
  filename.all <- paste(output.dir,"allnmprep_earlymdv_",sep="/")
  write.csv(orgdataLOQ, file=paste0(filename.org,"LOQ.csv"), quote=FALSE,row.names=FALSE)
  write.csv(orgdataHLOQ, file=paste0(filename.org,"HLOQ.csv"), quote=FALSE,row.names=FALSE)
  write.csv(alldataLOQ, file=paste0(filename.all,"LOQ.csv"), quote=FALSE,row.names=FALSE)
  write.csv(alldataHLOQ, file=paste0(filename.all,"HLOQ.csv"), quote=FALSE,row.names=FALSE)
