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

# -----------------------------------------------------------------------------
# Plot an Uppsala-style pcVPC
# -----------------------------------------------------------------------------
# PRED Correction
# Bergstrand et al 2011 - Prediction-Corrected Visual Predictive Checks for
# Diagnosing Nonlinear Mixed-Effects Model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate the median PRED for each TADBIN
  SIM.data1$PRED <- as.numeric(SIM.data1$PRED)
  SIM.data1BIN <- summaryBy(PRED ~ TADBIN, SIM.data1, FUN = median, na.rm = T)
  SIM.data2$PRED <- as.numeric(SIM.data2$PRED)
  SIM.data2BIN <- summaryBy(PRED ~ TADBIN, SIM.data2, FUN = median, na.rm = T)
	SIM.data3$PRED <- as.numeric(SIM.data3$PRED)
	SIM.data3BIN <- summaryBy(PRED ~ TADBIN, SIM.data3, FUN = median, na.rm = T)

# Merge median PREDs into simulated dataset matching for their TIMEBIN
	SIM.data1 <- merge(SIM.data1, SIM.data1BIN, by = c("TADBIN"), all = T)
  names(SIM.data1)[names(SIM.data1) == "PRED.median"] <- "PREDMED"
	SIM.data1 <- SIM.data1[with(SIM.data1,
		order(SIM.data1$SIM, SIM.data1$ID, SIM.data1$TAD, SIM.data1$TADBIN)), ]
	ORG.data1 <- ORG.data[with(ORG.data,
		order(ORG.data$ID, ORG.data$TAD, ORG.data$TADBIN)), ]

  SIM.data2 <- merge(SIM.data2, SIM.data2BIN, by = c("TADBIN"), all = T)
  names(SIM.data2)[names(SIM.data2) == "PRED.median"] <- "PREDMED"
	SIM.data2 <- SIM.data2[with(SIM.data2,
		order(SIM.data2$SIM, SIM.data2$ID, SIM.data2$TAD, SIM.data2$TADBIN)), ]
	ORG.data2 <- ORG.data[with(ORG.data,
		order(ORG.data$ID, ORG.data$TAD, ORG.data$TADBIN)), ]

  SIM.data3 <- merge(SIM.data3, SIM.data3BIN, by = c("TADBIN"), all = T)
  names(SIM.data3)[names(SIM.data3) == "PRED.median"] <- "PREDMED"
	SIM.data3 <- SIM.data3[with(SIM.data3,
		order(SIM.data3$SIM, SIM.data3$ID, SIM.data3$TAD, SIM.data3$TADBIN)), ]
	ORG.data3 <- ORG.data[with(ORG.data,
		order(ORG.data$ID, ORG.data$TAD, ORG.data$TADBIN)), ]

# Subset for one simulation of the same length of the original dataset
  SIM.data1ONE <- SIM.data1[SIM.data1$SIM == 1, ]
  SIM.data2ONE <- SIM.data2[SIM.data2$SIM == 1, ]
	SIM.data3ONE <- SIM.data3[SIM.data3$SIM == 1, ]

# Add median PRED for each TIMEBIN to the orignal dataset
	ORG.data1$PREDMED <- SIM.data1ONE$PREDMED
	ORG.data1$PRED <- SIM.data1ONE$PRED
  ORG.data2$PREDMED <- SIM.data2ONE$PREDMED
  ORG.data2$PRED <- SIM.data2ONE$PRED
  ORG.data3$PREDMED <- SIM.data3ONE$PREDMED
  ORG.data3$PRED <- SIM.data3ONE$PRED

# Calculate the prediction corrected observed and simulated DVs
	ORG.data1$pcY <- (ORG.data1$DV)*(ORG.data1$PREDMED)/(ORG.data1$PRED)
	SIM.data1$pcY <- (SIM.data1$DV)*(SIM.data1$PREDMED)/(SIM.data1$PRED)
  ORG.data2$pcY <- (ORG.data2$DV)*(ORG.data2$PREDMED)/(ORG.data2$PRED)
  SIM.data2$pcY <- (SIM.data2$DV)*(SIM.data2$PREDMED)/(SIM.data2$PRED)
  ORG.data3$pcY <- (ORG.data3$DV)*(ORG.data3$PREDMED)/(ORG.data3$PRED)
  SIM.data3$pcY <- (SIM.data3$DV)*(SIM.data3$PREDMED)/(SIM.data3$PRED)

# Combine the 3 SIM.data's and 3 ORG.datat's into a single SIM.data and ORG.data
  SIM.data1$MODEL <- 1
  SIM.data2$MODEL <- 2
  SIM.data3$MODEL <- 3
  SIM.data <- rbind(SIM.data1, SIM.data2, SIM.data3)
  SIM.data$MODEL <- factor(SIM.data$MODEL)
  levels(SIM.data$MODEL) <- c("Hughes et. al", "Guglieri-Lopez et. al", "Connarn et. al")

  ORG.data1$MODEL <- 1
  ORG.data2$MODEL <- 2
  ORG.data3$MODEL <- 3
  ORG.data <- rbind(ORG.data1, ORG.data2, ORG.data3)
  ORG.data$MODEL <- factor(ORG.data$MODEL)
  levels(ORG.data$MODEL) <- c("Hughes et. al", "Guglieri-Lopez et. al", "Connarn et. al")

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
		data.frame(
			medianS = median(x$pcY),
			loCI90S = CI90lo(x$pcY),
			hiCI90S = CI90hi(x$pcY)
		)
	})

# Build plot object
	titletext <- "VPC - Uppsala Style\n"
  p <- NULL
	p <- ggplot(data = ORG.data)
	# p <- p + ggtitle(titletext)

  p <- p + geom_point(aes(x = TADBIN, y = pcY, group = MODEL), colour = "blue", shape = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = pcY, group = MODEL), fun.y = median,
    geom = "line", colour = "red", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = pcY, group = MODEL), fun.y = CI90lo,
    geom = "line", colour = "red", linetype = "dashed", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = pcY, group = MODEL), fun.y = CI90hi,
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

	p <- p + scale_y_continuous("Prediction Corrected\nConcentration (mg/L)\n", lim = c(-0.05, 0.25))
	p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:8*3)
  p <- p + facet_wrap(~MODEL, ncol = 3)
	p

  ggsave("extval_pcvpc2.png", width = 24, height = 16, units = c("cm"))
