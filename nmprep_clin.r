###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()
   
# Set the working directory
  master.dir <- "D:/Hughes/Data"
  scriptname <- "nmprep_clin"
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
  source("D:/Hughes/functions_utility.r")
   
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
  file06003 <- "RAW_Clinical/datacheck_clin_06003_Output/06003_finaldata.csv"
	file05115 <- "RAW_Clinical/datacheck_clin_05115_Output/05115_finaldata.csv"
	file08056 <- "RAW_Clinical/datacheck_clin_08056_Output/08056_finaldata.csv"
	file10016 <- "RAW_Clinical/datacheck_clin_10016_Output/10016_finaldata.csv"
  data06003 <- read.csv(file06003, stringsAsFactors=F)
	data05115 <- read.csv(file05115, stringsAsFactors=F)
	data08056 <- read.csv(file08056, stringsAsFactors=F)
	data10016 <- read.csv(file10016, stringsAsFactors=F)
  
#Create nmprep data file
  datanew <- rbind(data06003,data05115,data08056,data10016)
	filename.out <- paste(output.dir,"fulldata.csv",sep="/")
  write.csv(datanew, file=filename.out, row.names=FALSE)
	
#Nmprep FINAL datacheck
	names(datanew)
  str(datanew)
  npat <- length(unique(datanew$X.ID))
  npat
	
#Week 1 only - Caucasian or Non-Caucasian
#ID TIME AMT EVID DV MDV AGE WT HT GEND RACE SECR DXCAT

  nmcols <- dataFIX[-c(2,3,4,6,8,10,13,18,19,21)]  #All columns except EVID
  nmcols$EVID <- 1
  nmcols$EVID[is.na(nmcols$AMT)] <- 0

  nmprep1 <- nmcols[datanew2$XSAMP==0,c(1,4,3,16,5,6,14,15,7,9,10,8,12,13,11,2)]
  nmprep1[is.na(nmprep1)] <- "."
  
  filename.out <- paste(output.dir,"06003_clin_nmprepwk1.csv",sep="/")
  write.csv(nmprep1, file=filename.out, quote=FALSE,row.names=FALSE)
  
#All Weeks
#ID TIME AMT EVID DV MDV ADDL II AGE WT HT GEND RACE SECR DXCAT
  nmprep2 <- nmcols[c(1,4,3,16,5,6,14,15,7,9,10,8,12,13,11,2)]
  nmprep2[is.na(nmprep2)] <- "."
  
  filename.out <- paste(output.dir,"06003_clin_nmprep.csv",sep="/")
  write.csv(nmprep1, file=filename.out, quote=FALSE,row.names=FALSE)
  
  simprep2 <- nmprep2
  simprep2$DV <- "."
  
  filename.out <- paste(output.dir,"06003_clin_simprep.csv",sep="/")
  write.csv(simprep2, file=filename.out, quote=FALSE,row.names=FALSE)