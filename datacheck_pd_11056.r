# Data check for PD data
# Aims are to determine what PD values are available then formulate a
# relationship between available values and SeCr which doesn't exist for this
# set of patients
# -----------------------------------------------------------------------------
# Set up environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_pd_11056"
  setwd(master.dir)

# Load libraries
  library(readxl)
  #library(vwr)  # visual word recognition
  library(plyr)
  library(ggplot2)

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
# -----------------------------------------------------------------------------
# Load in data
  file.name.in <- "RAW_Clinical/rawdata-lena_10156/Data_PK Analysis.xlsx"
  data.demo <- read_excel(file.name.in, sheet = 1)  # demographics
  # data.ig <- read_excel(file.name.in, sheet = 2)  # immunoglobulins
  data.tox <- read_excel(file.name.in, sheet = 3)  # toxicities
  data.res <- read_excel(file.name.in, sheet = 4)  # clinical response
  # data.sero <- read_excel(file.name.in, sheet = 5)  # serotype
# -----------------------------------------------------------------------------
# First clean up data.demo
  data.pd <- as.data.frame(data.demo)
  names(data.pd) <- c("SID", "ID", "AGE", "SEX", "WT", "HT")

# Now for the hard bits
# data.tox
  data.se <- as.data.frame(data.tox)
  names(data.se) <- c("SID", "CYC", "DAY", "TOX", "GRADE", "ONSET", "ATTR")

# Create ID column by splitting strings
  data.se.id <- as.numeric(
    unlist(
      strsplit(data.se$SID, "[^0-9]+")
    )  # unlist
  )  # as.numeric
  data.se$ID <- data.se.id[1:length(data.se$SID)*2]
  with(data.se, table(ID, useNA = "always"))

# Check the other columns
  with(data.se, table(CYC, useNA = "always"))
  with(data.se, table(DAY, useNA = "always"))
  with(data.se, table(TOX, useNA = "always"))
  # many categories have multiple entries due to spelling
  with(data.se, table(GRADE, useNA = "always"))
  with(data.se, table(ONSET, useNA = "always"))
  with(data.se, table(ATTR, useNA = "always"))
  # majority of toxicities are considered unrelated to lenalidomide
  # they are still important though, what attributes does not matter

# Firstly need to deal with the issue of inconsistent categories
  data.se.tox <- tolower(data.se$TOX)
  c(length(data.se$TOX), length(unique(data.se.tox)))
  # that was a massive improvement
  which(table(data.se.tox) > 5)
  # but there are still alot of categories with low numbers
  # and high numbers in categories which are effectively the same
  # oh wow, one of them is increased creatinine

# Toxicity seems to be mainly acute not pre-existing, so may not be ideal for 
# predicting SeCr...
# But concentrations are taken from the SECOND cycle for this trial..
# anything going on around that cycle might be useful, so lets go first 3 cycles
  se.sub <- data.se[data.se$CYC %in% 1:3, ]
  data.se.tox <- tolower(se.sub$TOX)
  c(length(se.sub$TOX), length(unique(data.se.tox)))

  data.se.tox[which(RecordLinkage::jarowinkler("increased creatinine", data.se.tox) > 0.8)]
  which(vwr::levenshtein.distance("creatinine", data.se.tox) < 7)
