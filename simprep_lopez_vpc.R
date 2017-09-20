# Prepare dataset for simulation
# -----------------------------------------------------------------------------
# Dataset will then be simulated using my model and the Guglieri-Lopez model
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/LOPEZ_SIM"
  setwd(master.dir)

# Load required packages
  library(ggplot2)
  library(Hmisc)

# Source functions
	source("E:/Hughes/functions_utility.r")

# Customize ggplot2 theme - R 2.15.3
  theme_bw2 <- theme_set(theme_bw(base_size = 22))
  theme_bw2 <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"),
    axis.title.x=element_text(size = 18, vjust = 0),
    axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
    strip.text.x=element_text(size = 16),
    strip.text.y=element_text(size = 16, angle = 90))

# Confidence intervals - from function utility
  CI90lo <- function(x) quantile(x, probs=0.05)
  CI90hi <- function(x) quantile(x, probs=0.95)

  CI95lo <- function(x) quantile(x, probs=0.025)
  CI95hi <- function(x) quantile(x, probs=0.975)

# -----------------------------------------------------------------------------
# Process the simulated *.fit file.
# This simulated file will be the "observed" data
  runname <- "RUN002_LOPEZ"
  # processSIMdata(paste(runname, ".ctl", sep = ""))
  ORG.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  ORG.data <- ORG.data[ORG.data$MDV == 0, ]

# This simulated file will be the simulation data
	runname <- "RUN001_BASE"
  # processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data <- SIM.data[SIM.data$MDV == 0, ]

# -----------------------------------------------------------------------------
# Assign factors to covariates
# Time binning
  timebin_cuts <- c(1.02, 2.02, 3.02, 4.02, 6.02, 8.02, 12.02, 16.02, 20.02, 24.02)
  ORG.data$TIMEBIN <- cut2(ORG.data$TIME, cuts = timebin_cuts, levels.mean = T)
  ORG.data$TIMEBIN <- as.numeric(paste(ORG.data$TIMEBIN))

  SIM.data$TIMEBIN <- cut2(SIM.data$TIME, cuts = timebin_cuts, levels.mean = T)
  SIM.data$TIMEBIN <- as.numeric(paste(SIM.data$TIMEBIN))

# Covariates
  ORG.data$IDf <- as.factor(ORG.data$ID)
  ORG.data$SEXf <- factor(ORG.data$SEX, labels=c("F", "M"))
  ORG.data$DOSEf <- factor(ifelse(ORG.data$DOSELVL<=5,1,2),
    labels=c("Dose <10mg","Dose >10mg"))
  ORG.data$CRCLf <- factor(ifelse(ORG.data$CRCL2<=60,1,2),
    labels=c("CrCl <60mL/min","CrCl >60mL/min"))

  SIM.data$IDf <- as.factor(SIM.data$ID)
  SIM.data$SEXf <- factor(SIM.data$SEX, labels=c("F", "M"))
  SIM.data$DOSEf <- factor(ifelse(SIM.data$DOSELVL<=5,1,2),
    labels=c("Dose <10mg","Dose >10mg"))
  SIM.data$CRCLf <- factor(ifelse(SIM.data$CRCL2<=60,1,2),
    labels=c("CrCl <60mL/min","CrCl >60mL/min"))

# -----------------------------------------------------------------------------
  # Uppsala Style VPC
  # Calculate 5, 50 and 95 percentiles for each simulated study (S)
  SIM.data.bystudy.median <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) median(df$DV))
  SIM.data.bystudy.median <- rename(SIM.data.bystudy.median, c("V1"="medianS"))

  SIM.data.bystudy.loCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90lo(df$DV))
  SIM.data.bystudy.loCI <- rename(SIM.data.bystudy.loCI, c("5%"="loCI90S"))

  SIM.data.bystudy.hiCI <- ddply(SIM.data, .(SIM,TIMEBIN), function(df) CI90hi(df$DV))
  SIM.data.bystudy.hiCI <- rename(SIM.data.bystudy.hiCI, c("95%"="hiCI90S"))

  SIM.data.bystudy <- data.frame(SIM.data.bystudy.median, "loCI90S"=SIM.data.bystudy.loCI$loCI90S, "hiCI90S"=SIM.data.bystudy.hiCI$hiCI90S)

  plotobj <- NULL
  plotobj2 <- NULL

  #Titletext with 3 lines
  titletext <- expression(atop(VPC~plot~Uppsala~Style,
                               atop(italic("Red solid line shows mean of data, dashed lines shows 90% CI of data"),
                                    "Black lines show median and 90% CI for simulated data. Ribbons show 95% CI around sim. data")))

  #titletext <- "VPC - Uppsala Style\n"
  plotobj <- ggplot(data=ORG.data)

  #Median simulated with confidence band
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="red")
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=medianS), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", size=1)

  #Lower 90% CI simulated with confidence band
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=loCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)

  #Upper 90% CI simulated with confidence band
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, geom="ribbon", fun.ymin="CI95lo", fun.ymax="CI95hi", alpha=0.3, fill="blue")
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=hiCI90S), data=SIM.data.bystudy, fun.y=median, geom="line", colour="black", linetype="dashed", size=1)


  plotobj <- plotobj + geom_point(aes(x=TIMEBIN, y=DV), colour="blue", shape = 1)
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=median, geom="line", colour="red", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90lo, geom="line", colour="red", linetype="dashed", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), fun.y=CI90hi, geom="line", colour="red", linetype="dashed", size=1)
  # For simulated data: add median, CI90lo, CI90hi
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=median, geom="line", colour="black", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90lo, geom="line", colour="black", size=1)
  plotobj <- plotobj + stat_summary(aes(x=TIMEBIN, y=DV), data=SIM.data, fun.y=CI90hi, geom="line", colour="black", size=1)

  plotobj <- plotobj +  theme(plot.title = element_text(size = rel(1)))
  plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + scale_y_continuous("Concentration (ug/L)\n")
  plotobj <- plotobj + scale_x_log10("\nTime after dose (hours)")
  #plotobj <- plotobj + facet_wrap(~DVIDf)
  plotobj <- plotobj + theme(strip.background = element_rect(fill = "grey95", colour = "grey50"))
  plotobj

  #normal scale- non-facetted
  ggsave("Uppsala_VPC.png", width=20, height=16, units=c("cm"))

  # log Y-scale-non- facetted
  plotobj2 <- plotobj + scale_y_log10("Concentration (ug/L)\n")
  plotobj2
  ggsave("Uppsala_VPC_log.png", width=20, height=16, units=c("cm"))
