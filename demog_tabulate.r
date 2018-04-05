# Tabulate the demographic information from the studies for publishing
# -----------------------------------------------------------------------------
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "demog_tabulate"
  setwd(master.dir)

# Load libraries
  library(plyr)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

# Organise working and output directories
  working.dir <- paste(master.dir,"RAW_Clinical",sep="/")
  workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

  output.dir <- paste(working.dir,"/",scriptname,"_Output",sep="")
  if(!file.exists(output.dir)){
	 dir.create(output.dir)
  }

# -----------------------------------------------------------------------------
# Load the data
  study.numbers <- c("05115", "06003", "08056", "10016", "10156")
  file.dir <- paste0(working.dir, "/datacheck_clin_", study.numbers, "_Output")
  file.names <- paste0(file.dir, "/", study.numbers, "_covdata.csv")

  datain <- llply(file.names, function(x) {
    read.csv(x, stringsAsFactors = F, na.strings = c("."))
  })

# -----------------------------------------------------------------------------
# Calculate summary statistics
  data <- llply(datain, function(x) {
    if (length(x$SECRMGDL) == 0) {
      x$CRCL <- 0
    } else {
      fsex <- ifelse(x$GEND==1,1.23,1.04)
      ibw <- ifelse(x$GEND==1,50+0.9*(x$HT-152),45.5+0.9*(x$HT-152))
      x$CRCL <- (140-x$AGE)*ibw*fsex/(x$SECRMGDL*88.4)
    }
    ffm_c1 <- x$GEND*6.68e3
    ffm_c1[ffm_c1 == 0] <- 8.78e3
    ffm_c2 <- x$GEND*216
    ffm_c2[ffm_c2 == 0] <- 244
    bmi <- x$WT/(x$HT/100)^2
    x$FFM <- 9.27e3 * x$WT / (ffm_c1 + ffm_c2 * bmi)
    x
  })

  ldply(data, function(x) {
    c(no = length(x$ID), nsex = length(x$GEND[x$GEND == 0]), msex = 1 - mean(x$GEND),
      meanage = mean(x$AGE), sdage = sd(x$AGE),
      minage = min(x$AGE), medage = median(x$AGE), maxage = max(x$AGE),
      meanwt = mean(x$WT), sdwt = sd(x$WT),
      minwt = min(x$WT), medwt = median(x$WT), maxwt = max(x$WT),
      meanwt = mean(x$HT), sdht = sd(x$WT),
      minht = min(x$HT), medht = median(x$HT), maxht = max(x$HT),
      meansecr = mean(x$SECRMGDL*88.4), sdsecr = sd(x$SECRMGDL*88.4),
      minsecr = min(x$SECRMGDL*88.4), medsecr = median(x$SECRMGDL*88.4), maxsecr = max(x$SECRMGDL*88.4),
      meancrcl = mean(x$CRCL), sdcrcl = sd(x$CRCL),
      mincrcl = min(x$CRCL), medcrcl = median(x$CRCL), maxcrcl = max(x$CRCL),
      meanffm = mean(x$FFM), sdffm = sd(x$FFM),
      minffm = min(x$FFM), medffm = median(x$FFM), maxffm = max(x$FFM),
      meandose = mean(x$DOSEMG), sddose = sd(x$DOSEMG),
      mindose = min(x$DOSEMG), meddose = median(x$DOSEMG), maxdose = max(x$DOSEMG)
    )
  })

# -----------------------------------------------------------------------------
  all.studies <- rbind(data[[1]], data[[2]], data[[3]], data[[4]])
  sum.fun.all <- function(x) {
    c(no = length(x$ID), nsex = length(x$GEND[x$GEND == 0]), msex = 1 - mean(x$GEND),
      meanage = mean(x$AGE), sdage = sd(x$AGE),
      minage = min(x$AGE), medage = median(x$AGE), maxage = max(x$AGE),
      meanwt = mean(x$WT), sdwt = sd(x$WT),
      minwt = min(x$WT), medwt = median(x$WT), maxwt = max(x$WT),
      meanwt = mean(x$HT), sdht = sd(x$WT),
      minht = min(x$HT), medht = median(x$HT), maxht = max(x$HT),
      meansecr = mean(x$SECRMGDL*88.4), sdsecr = sd(x$SECRMGDL*88.4),
      minsecr = min(x$SECRMGDL*88.4), medsecr = median(x$SECRMGDL*88.4), maxsecr = max(x$SECRMGDL*88.4),
      meancrcl = mean(x$CRCL), sdcrcl = sd(x$CRCL),
      mincrcl = min(x$CRCL), medcrcl = median(x$CRCL), maxcrcl = max(x$CRCL),
      meanffm = mean(x$FFM), sdffm = sd(x$FFM),
      minffm = min(x$FFM), medffm = median(x$FFM), maxffm = max(x$FFM),
      meandose = mean(x$DOSEMG), sddose = sd(x$DOSEMG),
      mindose = min(x$DOSEMG), meddose = median(x$DOSEMG), maxdose = max(x$DOSEMG)
    )
  }
  sum.fun.all(all.studies)
  with(all.studies, table(DXCATNUM))
  with(all.studies, table(DOSEMG))

  maledata <- ldply(data[1:4], function(x) {
    indata <- x[x$GEND == 1,]
    out <- indata[,c("UID", "STUDY", "AGE", "WT", "HT")]
  })

  femaledata <- ldply(data[1:4], function(x) {
    indata <- x[x$GEND == 0,]
    out <- indata[,c("UID", "STUDY", "AGE", "WT", "HT")]
  })
