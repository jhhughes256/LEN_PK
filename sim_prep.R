# sim_prep.r

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "sim_prep"
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

# ------------------------------------------------------------------------------
  times <- unique(c(seq(0, 4, 0.25), seq(4, 8, 0.5), seq(8, 24, 1)))
  nobs <- length(times)

  simprep <- data.frame(
    ID = 1,
    TIME = times,
    AMT = c(25, rep(".", nobs - 1)),
    EVID = c(1, rep(0, nobs - 1)),
    DV = c(".", rep(1, nobs - 1)),
    CMT = c(1, rep(2, nobs - 1)),
    MDV = c(1, rep(0, nobs - 1)),
    FFM = 55
  )

  names(simprep)[1] <- "#ID"
  filename.out <- paste(output.dir,"simdata.csv",sep="/")
  write.csv(simprep, file = filename.out, quote = F, row.names = F)

# ------------------------------------------------------------------------------
#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()


#Set the working directory - to the parent directory where you do the modelling/VPC
	master.dir<-"E:/Hughes/Data/RAW_Clinical/sim_prep_Output"
	setwd(master.dir)

#Load libraries
	library(R2HTML)
	library(ggplot2)
	library(doBy)
	library(stringr)
	library(Hmisc)
	library(grid)
	library(plyr)

#Source WfNviaR functions file - Need to have this loaded as it runs most of the processing
	source("E:/Hughes/functions_utility.r")

#Use custom ggplot 2 theme
	theme_custom <- theme_set(theme_bw(18))
	theme_custom <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"))

#-------------------------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
  theme_bw2 <- theme_set(theme_bw(base_size = 22))
  theme_bw2 <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"),
                           axis.title.x=element_text(size = 18, vjust = 0),
                           axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
                           strip.text.x=element_text(size = 16),
                           strip.text.y=element_text(size = 16, angle = 90))
#-------------------------------------------------------------------------------
#Confidence intervals - from function utility
  CI90lo <- function(x) quantile(x, probs=0.05)
  CI90hi <- function(x) quantile(x, probs=0.95)

  CI95lo <- function(x) quantile(x, probs=0.025)
  CI95hi <- function(x) quantile(x, probs=0.975)

#-------------------------------------------------------------------------------
#Process the simulated *.fit file.
#Run name - Change this to the RUN you want to process
	runname <- "sims_1000"

#Process the fit file - Comment this out if you have already generated the csv; this will save time!
  processSIMdata(paste(runname,".ctl",sep=""))    # from the FUNCTION UTILITY

#Read the simulated data
  SIM.data <- read.csv(paste(runname,".nm7/",runname,".fit.csv",sep=""), stringsAsFactors=F, na.strings=".")

  plotobj <- ggplot(data=SIM.data)

  plotobj <- plotobj + stat_summary(aes(x=TIME, y=IPRED), fun.ymin=CI90lo, fun.ymax=CI90hi, geom="ribbon", fill="blue", alpha = 0.3)
  plotobj <- plotobj + stat_summary(aes(x=TIME, y=IPRED), fun.y=median, geom="line", colour="red", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIME, y=IPRED), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIME, y=IPRED), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)

  plotobj <- plotobj + theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)\n")
  plotobj <- plotobj + scale_x_log10("\nTime after dose (hours)")
  #plotobj <- plotobj + facet_wrap(~DVIDf)
  plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj
