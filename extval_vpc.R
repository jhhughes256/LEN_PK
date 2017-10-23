# Generating VPCs for each of the three models using the external validation dataset
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/EXT_VAL"
  setwd(master.dir)

# Load required packages
  library(ggplot2)
  library(Hmisc)
	library(doBy)
	library(plyr)

# Source functions
	source("E:/Hughes/functions_utility.r")

# Customize ggplot2 theme - R 2.15.3
  theme_bw2 <- theme_set(theme_bw(base_size = 22))
  theme_bw2 <- theme_update(plot.margin = unit(c(1, 0.5, 3, 0.5), "lines"),
    axis.title.x = element_text(size = 18, vjust = 0),
    axis.title.y = element_text(size = 18, vjust = 0, angle = 90),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16, angle = 90))

# Confidence intervals - from function utility
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)

  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)

# -----------------------------------------------------------------------------
# Read in data for plotting
# Process the simulated *.fit files
	runname1 <- "RUN001_BASE"
  # processSIMdata(paste(runname1,".ctl",sep=""))
  SIM.data1 <- read.csv(paste(runname1, ".nm7/", runname1, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data1 <- SIM.data1[SIM.data1$MDV == 0, ]

  runname2 <- "RUN002_LOPEZ"
  # processSIMdata(paste(runname2,".ctl",sep=""))
  SIM.data2 <- read.csv(paste(runname2, ".nm7/", runname2, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data2 <- SIM.data2[SIM.data2$MDV == 0, ]

  runname3 <- "RUN003_CELGENE"
  # processSIMdata(paste(runname3,".ctl",sep=""))
  SIM.data3 <- read.csv(paste(runname3, ".nm7/", runname3, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data3 <- SIM.data3[SIM.data3$MDV == 0, ]
  SIM.data3$DV <- exp(SIM.data3$DV)

# Read in the original data
  ORG.data <- read.csv("extval_10156.csv", stringsAsFactors = F, na.strings = ".")
	names(ORG.data)[names(ORG.data) == "X.ID"] <- "ID"
  ORG.data <- ORG.data[ORG.data$MDV == 0, ]

# -----------------------------------------------------------------------------
# Assign factors to covariates
# Time binning
  bin_cuts <- c(0.25, 0.75, 1.5, 2.5, 3.75, 5.5, 7, 25)
  ORG.data$TADBIN <- cut2(ORG.data$TAD, cuts = bin_cuts, levels.mean = T)
  ORG.data$TADBIN <- as.numeric(paste(ORG.data$TADBIN))
  # with(ORG.data, table(TADBIN))

  SIM.data1$TADBIN <- cut2(SIM.data1$TAD, cuts = bin_cuts, levels.mean = T)
  SIM.data1$TADBIN <- as.numeric(paste(SIM.data1$TADBIN))
  SIM.data2$TADBIN <- cut2(SIM.data2$TAD, cuts = bin_cuts, levels.mean = T)
  SIM.data2$TADBIN <- as.numeric(paste(SIM.data2$TADBIN))
  SIM.data3$TADBIN <- cut2(SIM.data3$TAD, cuts = bin_cuts, levels.mean = T)
  SIM.data3$TADBIN <- as.numeric(paste(SIM.data3$TADBIN))

# Covariates
	ORG.data$IDf <- as.factor(ORG.data$ID)
	ORG.data$SEXf <- factor(ORG.data$SEX, labels = c("F", "M"))
	ORG.data$CRCLf <- factor(ifelse(ORG.data$CRCL2 <= 60, 1, 2),
	  labels = c("CrCl <60mL/min", "CrCl >60mL/min"))

  SIM.data1$IDf <- as.factor(SIM.data1$ID)
  SIM.data2$IDf <- as.factor(SIM.data2$ID)
	SIM.data3$IDf <- as.factor(SIM.data3$ID)
  SIM.data1$SEXf <- factor(SIM.data1$SEX, labels = c("F", "M"))
  SIM.data2$SEXf <- factor(SIM.data2$SEX, labels = c("F", "M"))
	SIM.data3$SEXf <- factor(SIM.data3$SEX, labels = c("F", "M"))
  SIM.data1$CRCLf <- factor(ifelse(SIM.data1$CRCL2 <= 60, 1, 2),
    labels = c("CrCl <60mL/min", "CrCl >60mL/min"))
  SIM.data2$CRCLf <- factor(ifelse(SIM.data2$CRCL2 <= 60, 1, 2),
    labels = c("CrCl <60mL/min", "CrCl >60mL/min"))
	SIM.data3$CRCLf <- factor(ifelse(SIM.data3$CRCL2 <= 60, 1, 2),
	  labels = c("CrCl <60mL/min", "CrCl >60mL/min"))

# Combine the 3 SIM.data's into one!
  SIM.data1$MODEL <- 1
  SIM.data2$MODEL <- 2
  SIM.data3$MODEL <- 3
  SIM.data <- rbind(SIM.data1, SIM.data2, SIM.data3)
  SIM.data$MODEL <- factor(SIM.data$MODEL)
  levels(SIM.data$MODEL) <- c("Hughes et. al", "Guglieri-Lopez et. al", "Connarn et. al")

# -----------------------------------------------------------------------------
# Create pcVPC using Xpose method
# Plot the confidence interval for the simulated data's percentiles for each bin
# (for each simulated data set compute the percentiles for each bin, then, from
# all of the percentiles from all of the simulated datasets compute the 95% CI
# of these percentiles).
# http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate 5, 50 and 95 percentiles for each simulated study (S)
  SIM.data.bystudy <- ddply(SIM.data, .(MODEL, SIM, TADBIN), function(x) {
    dv <- x$DV
    data.frame(
      medianS = median(dv),
      loCI90S = CI90lo(dv),
      hiCI90S = CI90hi(dv)
    )
  })

# Build plot object
	titletext <- "VPC - Uppsala Style\n"
  p <- NULL
	p <- ggplot(data = ORG.data)
	p <- p + ggtitle(titletext)

  p <- p + geom_point(aes(x = TADBIN, y = DV), colour = "blue", shape = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = DV), fun.y = median,
    geom = "line", colour = "red", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = DV), fun.y = CI90lo,
    geom = "line", colour = "red", linetype = "dashed", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = DV), fun.y = CI90hi,
    geom = "line", colour = "red", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = medianS, group = MODEL), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
  p <- p + stat_summary(aes(x = TADBIN, y = medianS, group = MODEL), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S, group = MODEL), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S, group = MODEL), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S, group = MODEL), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S, group = MODEL), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + scale_y_continuous("Concentration (mg/L)\n")
	p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:12*2)
  p <- p + facet_wrap(~MODEL, ncol = 3)
	p
