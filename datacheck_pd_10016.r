###datacheck.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 10th July 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "datacheck_pd_10016"
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
#                              Prev. Treated BM                                #
  file.name.in <- "RAW_Clinical/rawdata-lena_10016_prev_treated_BM.xls"
  data.all <- read_excel(file.name.in, sheet=1)  #multiple dataframes within one sheet

  ct.sum.ptbm.raw <- data.all[2:73,11:20]
  names(ct.sum.ptbm.raw) <- c("SID","PKID","Day","CTval","GAPDH","Pgp","BCRP","CEBPa","CRBN")
  ct.sum.ptbm <- data.frame(
    "SID" = repeat.before(ct.sum.ptbm.raw$SID),
    "PKID" = as.numeric(repeat.before(ct.sum.ptbm.raw$PKID)),
    "Day" = as.numeric(repeat.before(str_extract(ct.sum.ptbm.raw$Day,"[0-9]"))),
    "CTval" = as.numeric(str_extract(ct.sum.ptbm.raw$CTval,"[0-9]")),
    "GAPDH" = as.numeric(ct.sum.ptbm.raw$GAPDH),
    "Pgp" = as.numeric(ct.sum.ptbm.raw$Pgp),
    "BCRP" = as.numeric(ct.sum.ptbm.raw$BCRP),
    "CEBPa" = as.numeric(ct.sum.ptbm.raw$CEBPa),
    "CRBN" = as.numeric(ct.sum.ptbm.raw$CRBN))

  dct.sum.ptbm <- ddply(ct.sum.ptbm, .(PKID, Day), function(x) {
    out <- data.frame(
      TRT = "PT",
      Tissue = "BM",
      Gene = c("Pgp","BCRP","CEBPa","CRBN"),
      GAPDHmean = rep(mean(x$GAPDH),4),
      GAPDHsd = rep(sd(x$GAPDH),4),
      DVmean = c(mean(x$Pgp), mean(x$BCRP), mean(x$CEBPa), mean(x$CRBN)),
      DVsd = c(sd(x$Pgp), sd(x$BCRP), sd(x$CEBPa), sd(x$CRBN)))

    out$dCt <- out$DVmean - out$GAPDHmean
    out$dCtmsd <- ((out$DVsd - out$GAPDHsd)^2)^0.5
    out
  })

  ddct.sum.ptbm <- ddply(dct.sum.ptbm, .(Gene), function(x) {
    max <- max(x$dCt, na.rm=TRUE)
    max.msd <- x$dCtmsd[which(x$dCt == max(x$dCt, na.rm=TRUE))]

    out <- x
    out$ddCt <- out$dCt - max
    out$ddCtmsd <- ((out$dCtmsd - max.msd)^2)^0.5
    out
  })

#                                Untreated BM                                  #
  file.name.in <- "RAW_Clinical/rawdata-lena_10016_untreated_BM.xls"
  data.all <- read_excel(file.name.in, sheet=1)  #multiple dataframes within one sheet

  ct.sum.utbm.raw <- data.all[2:45,11:20]
  names(ct.sum.utbm.raw) <- c("SID","PKID","Day","CTval","GAPDH","Pgp","BCRP","CEBPa","CRBN")
  ct.sum.utbm <- data.frame(
    "SID" = repeat.before(ct.sum.utbm.raw$SID),
    "PKID" = as.numeric(repeat.before(ct.sum.utbm.raw$PKID)),
    "Day" = as.numeric(repeat.before(str_extract(ct.sum.utbm.raw$Day,"[0-9]"))),
    "CTval" = as.numeric(str_extract(ct.sum.utbm.raw$CTval,"[0-9]")),
    "GAPDH" = as.numeric(ct.sum.utbm.raw$GAPDH),
    "Pgp" = as.numeric(ct.sum.utbm.raw$Pgp),
    "BCRP" = as.numeric(ct.sum.utbm.raw$BCRP),
    "CEBPa" = as.numeric(ct.sum.utbm.raw$CEBPa),
    "CRBN" = as.numeric(ct.sum.utbm.raw$CRBN))

  dct.sum.utbm <- ddply(ct.sum.utbm, .(PKID, Day), function(x) {
    out <- data.frame(
      TRT = "UT",
      Tissue = "BM",
      Gene = c("Pgp","BCRP","CEBPa","CRBN"),
      GAPDHmean = rep(mean(x$GAPDH),4),
      GAPDHsd = rep(sd(x$GAPDH),4),
      DVmean = c(mean(x$Pgp), mean(x$BCRP), mean(x$CEBPa), mean(x$CRBN)),
      DVsd = c(sd(x$Pgp), sd(x$BCRP), sd(x$CEBPa), sd(x$CRBN)))

    out$dCt <- out$DVmean - out$GAPDHmean
    out$dCtmsd <- ((out$DVsd - out$GAPDHsd)^2)^0.5
    out
  })

  ddct.sum.utbm <- ddply(dct.sum.utbm, .(Gene), function(x) {
    max <- max(x$dCt, na.rm=TRUE)
    max.msd <- x$dCtmsd[which(x$dCt == max(x$dCt, na.rm=TRUE))]

    out <- x
    out$ddCt <- out$dCt - max
    out$ddCtmsd <- ((out$dCtmsd - max.msd)^2)^0.5
    out
  })

#                            Prev. Treated PBMC                                #
file.name.in <- "RAW_Clinical/rawdata-lena_10016_prev_treated_PBMC.xls"
data.all <- read_excel(file.name.in, sheet=1)  #multiple dataframes within one sheet

ct.sum.ptpbmc.raw <- data.all[2:57,11:20]
names(ct.sum.ptpbmc.raw) <- c("SID","PKID","Day","CTval","GAPDH","Pgp","BCRP","CEBPa","CRBN")
ct.sum.ptpbmc <- data.frame(
  "SID" = repeat.before(ct.sum.ptpbmc.raw$SID),
  "PKID" = as.numeric(repeat.before(ct.sum.ptpbmc.raw$PKID)),
  "Day" = as.numeric(repeat.before(str_extract(ct.sum.ptpbmc.raw$Day,"[0-9]"))),
  "CTval" = as.numeric(str_extract(ct.sum.ptpbmc.raw$CTval,"[0-9]")),
  "GAPDH" = as.numeric(ct.sum.ptpbmc.raw$GAPDH),
  "Pgp" = as.numeric(ct.sum.ptpbmc.raw$Pgp),
  "BCRP" = as.numeric(ct.sum.ptpbmc.raw$BCRP),
  "CEBPa" = as.numeric(ct.sum.ptpbmc.raw$CEBPa),
  "CRBN" = as.numeric(ct.sum.ptpbmc.raw$CRBN))

dct.sum.ptpbmc <- ddply(ct.sum.ptpbmc, .(PKID, Day), function(x) {
  out <- data.frame(
    TRT = "PT",
    Tissue = "PBMC",
    Gene = c("Pgp","BCRP","CEBPa","CRBN"),
    GAPDHmean = rep(mean(x$GAPDH),4),
    GAPDHsd = rep(sd(x$GAPDH),4),
    DVmean = c(mean(x$Pgp), mean(x$BCRP), mean(x$CEBPa), mean(x$CRBN)),
    DVsd = c(sd(x$Pgp), sd(x$BCRP), sd(x$CEBPa), sd(x$CRBN)))

  out$dCt <- out$DVmean - out$GAPDHmean
  out$dCtmsd <- ((out$DVsd - out$GAPDHsd)^2)^0.5
  out
})

ddct.sum.ptpbmc <- ddply(dct.sum.ptpbmc, .(Gene), function(x) {
  max <- max(x$dCt, na.rm=TRUE)
  max.msd <- x$dCtmsd[which(x$dCt == max(x$dCt, na.rm=TRUE))]

  out <- x
  out$ddCt <- out$dCt - max
  out$ddCtmsd <- ((out$dCtmsd - max.msd)^2)^0.5
  out
})

#                              Untreated PBMC                                  #
file.name.in <- "RAW_Clinical/rawdata-lena_10016_untreated_PBMC.xls"
data.all <- read_excel(file.name.in, sheet=2)  #multiple dataframes within one sheet

ct.sum.utpbmc.raw <- data.all[3:58,11:20]
names(ct.sum.utpbmc.raw) <- c("SID","PKID","Day","CTval","GAPDH","Pgp","BCRP","CEBPa","CRBN")
ct.sum.utpbmc <- data.frame(
  "SID" = repeat.before(ct.sum.utpbmc.raw$SID),
  "PKID" = as.numeric(repeat.before(ct.sum.utpbmc.raw$PKID)),
  "Day" = as.numeric(repeat.before(str_extract(ct.sum.utpbmc.raw$Day,"[0-9]"))),
  "CTval" = as.numeric(str_extract(ct.sum.utpbmc.raw$CTval,"[0-9]")),
  "GAPDH" = as.numeric(ct.sum.utpbmc.raw$GAPDH),
  "Pgp" = as.numeric(ct.sum.utpbmc.raw$Pgp),
  "BCRP" = as.numeric(ct.sum.utpbmc.raw$BCRP),
  "CEBPa" = as.numeric(ct.sum.utpbmc.raw$CEBPa),
  "CRBN" = as.numeric(ct.sum.utpbmc.raw$CRBN))

dct.sum.utpbmc <- ddply(ct.sum.utpbmc, .(PKID, Day), function(x) {
  out <- data.frame(
    TRT = "UT",
    Tissue = "PBMC",
    Gene = c("Pgp","BCRP","CEBPa","CRBN"),
    GAPDHmean = rep(mean(x$GAPDH),4),
    GAPDHsd = rep(sd(x$GAPDH),4),
    DVmean = c(mean(x$Pgp), mean(x$BCRP), mean(x$CEBPa), mean(x$CRBN)),
    DVsd = c(sd(x$Pgp), sd(x$BCRP), sd(x$CEBPa), sd(x$CRBN)))

  out$dCt <- out$DVmean - out$GAPDHmean
  out$dCtmsd <- ((out$DVsd - out$GAPDHsd)^2)^0.5
  out
})

ddct.sum.utpbmc <- ddply(dct.sum.utpbmc, .(Gene), function(x) {
  max <- max(x$dCt, na.rm=TRUE)
  max.msd <- x$dCtmsd[which(x$dCt == max(x$dCt, na.rm=TRUE))]

  out <- x
  out$ddCt <- out$dCt - max
  out$ddCtmsd <- ((out$dCtmsd - max.msd)^2)^0.5
  out
})

### Combine data and save to file
  datapd <- rbind(ddct.sum.ptbm, ddct.sum.utbm, ddct.sum.ptpbmc, ddct.sum.utpbmc)
  filename.out <- paste(output.dir,"10016_allpd.csv",sep="/")
  write.csv(datapd, filename.out)

### Plot data
