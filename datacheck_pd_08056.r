###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_pd_08056"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)
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

### ------------------------------ PD Data --------------------------------- ###
  ## missingdata.csv
  file.name.in2 <- "RAW_Clinical/rawdata-lena_08056_missingdata.csv"
  datamis <- read.csv(file.name.in2,stringsAsFactors=F, na.strings=c("."))[(1:25)*2,1:2]  #take every second value
  names(datamis) <- c("SEQ","ID")

  ## demog_excel.xls - 13 sheets of data
  file.name.in3 <- "RAW_Clinical/rawdata-lena_08056_demog_excel.xlsx"
  #rows removed for double values; empty first column removed
  #demog.cat <- read_excel(file.name.in3, sheet=1)[-c(5,39),-1]  #demographxics
  #demog.offstudy <- read_excel(file.name.in3, sheet=2)  #removed from study?
  #demog.ps <- read_excel(file.name.in3, sheet=3)  #past surgery?
  #demog.pt <- read_excel(file.name.in3, sheet=4)  #past treatment
  #demog.bl <- read_excel(file.name.in3, sheet=5)  #toxicity around symptoms
  #demog.cont <- read_excel(file.name.in3, sheet=6)  #weight, height, bsa, with accompanying cycle and visit dates
  #demog.ca <- read_excel(file.name.in3, sheet=7)  #treatment information
  #demog.tox <- read_excel(file.name.in3, sheet=8)  #toxicity around bloods
  demog.haem <- read_excel(file.name.in3, sheet=9)  #haematology results - extensive
  #demog.secr <- read_excel(file.name.in3, sheet=10)  #creatinine levels - extensive
  #demog.ig <- read_excel(file.name.in3, sheet=11)  #ig levels
  #demog.ul <- read_excel(file.name.in3, sheet=12)  #upper lymph nodes?
  #demog.ll <- read_excel(file.name.in3, sheet=13)  #lower lymph nodes?

  names(demog.haem) <- c("SEQ","INIT","CYCLE","DAY","DATE","TIME","EXAM","DV")
  HCT.raw <- demog.haem[demog.haem$EXAM == "Hematocrit", ]
  HCT.all <- merge(datamis,HCT.raw)
  HCT.sub <- na.omit(HCT.all[HCT.all$CYCLE == 2 & as.numeric(HCT.all$DAY) < 10, ])

  HCT.mean <- ldply(seq_len(length(unique(HCT.sub$ID))), function(x) {
    sub <- HCT.sub[HCT.sub$ID == unique(HCT.sub$ID)[x], ]
    out <- data.frame(
    SEQ = sub$SEQ[1],
    ID = sub$ID[1],
    MEAN = mean(as.numeric(sub$DV)))
    out
  })

  write.csv(HCT.mean,paste(output.dir,"HCTmean.csv",sep="/"))
