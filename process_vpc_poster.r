# Generating VPCs for both the ACoP8 poster and publication
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/FLAG/COV15"
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
# Process the simulated *.fit file.
	runname <- "RUN038_CL_CRCL2_FFM_CATGT150_2_VPC"
  processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data <- SIM.data[SIM.data$MDV == 0, ]

# Read in the original data
  ORG.data <- read.csv("allnmprep_flagged.csv", stringsAsFactors = F, na.strings = ".")
	names(ORG.data)[names(ORG.data) == "X.ID"] <- "ID"
  ORG.data <- ORG.data[ORG.data$MDV == 0 & ORG.data$FLAG == 0, ]

# -----------------------------------------------------------------------------
# Assign factors to covariates
# Time binning
	bin_cuts <- c(0.52, 1.02, 2.02, 3.02, 5.02, 8.02, 49)
  ORG.data$TADBIN <- cut2(ORG.data$TAD, cuts = bin_cuts, levels.mean = T)
  ORG.data$TADBIN <- as.numeric(paste(ORG.data$TADBIN))
	# with(ORG.data, table(TADBIN))

  SIM.data$TADBIN <- cut2(SIM.data$TAD, cuts = bin_cuts, levels.mean = T)
  SIM.data$TADBIN <- as.numeric(paste(SIM.data$TADBIN))

# Covariates
	ORG.data$IDf <- as.factor(ORG.data$ID)
	ORG.data$RACEf <- factor(ORG.data$RACE, labels = c("White", "Other"))
	ORG.data$SEXf <- factor(ORG.data$SEX, labels = c("F", "M"))
	ORG.data$DXCATf <- factor(ORG.data$DXCATNUM, labels = c("CLL", "AML", "ALL", "MM"))
	ORG.data$STUDYf <- factor(ORG.data$STUDY)
	ORG.data$DOSEf <- factor(ifelse(ORG.data$DOSELVL <= 5, 1, 2),
	  labels = c("Dose <10mg", "Dose >10mg"))
	ORG.data$CRCLf <- factor(ifelse(ORG.data$CRCL2 <= 60, 1, 2),
	  labels = c("CrCl <60mL/min", "CrCl >60mL/min"))
  ORG.data$CRCL2f <- 1
  ORG.data$CRCL2f[ORG.data$CRCL2 >= 60] <- 2
  ORG.data$CRCL2f[ORG.data$CRCL2 >= 90] <- 3
  ORG.data$CRCL2f[ORG.data$CRCL2 > 120] <- 4
  ORG.data$CRCL2f <- factor(ORG.data$CRCL2f,
    labels = c("CrCl <60mL/min", "CrCl 60-90mL/min", "CrCl 90-120mL/min", "CrCl >120mL/min"))

	SIM.data$IDf <- as.factor(SIM.data$ID)
	SIM.data$RACEf <- factor(SIM.data$RACE, labels = c("White", "Other"))
	SIM.data$SEXf <- factor(SIM.data$SEX, labels = c("F", "M"))
	SIM.data$DXCATf <- factor(SIM.data$DXCAT, labels = c("CLL", "AML", "ALL", "MM"))
	SIM.data$STUDYf <- factor(SIM.data$STUDY)
	SIM.data$DOSEf <- factor(ifelse(SIM.data$DOSELVL <= 5, 1, 2),
	  labels = c("Dose <10mg", "Dose >10mg"))
	SIM.data$CRCLf <- factor(ifelse(SIM.data$CRCL2 <= 60, 1, 2),
	  labels = c("CrCl <60mL/min", "CrCl >60mL/min"))
  SIM.data$CRCL2f <- 1
  SIM.data$CRCL2f[SIM.data$CRCL2 >= 60] <- 2
  SIM.data$CRCL2f[SIM.data$CRCL2 >= 90] <- 3
  SIM.data$CRCL2f[SIM.data$CRCL2 > 120] <- 4
  SIM.data$CRCL2f <- factor(SIM.data$CRCL2f,
    labels = c("CrCl <60mL/min", "CrCl 60-90mL/min", "CrCl 90-120mL/min", "CrCl >120mL/min"))

# -----------------------------------------------------------------------------
# Create pcVPC using Xpose method
# Plot the confidence interval for the simulated data's percentiles for each bin
# (for each simulated data set compute the percentiles for each bin, then, from
# all of the percentiles from all of the simulated datasets compute the 95% CI
# of these percentiles).
# http://www.inside-r.org/packages/cran/xpose4specific/docs/xpose.VPC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate 5, 50 and 95 percentiles for each simulated study (S)
	SIM.data.bystudy <- ddply(SIM.data, .(SIM, TADBIN), function(x) {
		data.frame(
			medianS = median(x$DV),
			loCI90S = CI90lo(x$DV),
			hiCI90S = CI90hi(x$DV)
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

	p <- p + stat_summary(aes(x = TADBIN, y = medianS), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
  p <- p + stat_summary(aes(x = TADBIN, y = medianS), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S), data = SIM.data.bystudy,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S), data = SIM.data.bystudy,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + scale_y_log10("Concentration (mg/L)\n")
	p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:12*2)
	p

  ggsave(paste0(runname, ".nm7/VPC_poster.png"),
	  width = 20, height = 16, units = c("cm"))

  p + facet_wrap(~CRCL2f)
  ggsave(paste0(runname, ".nm7/VPC_byCRCL.png"),
    width = 20, height = 16, units = c("cm"))
