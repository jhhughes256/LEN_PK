# Generating VPCs for both the ACoP8 poster and publication
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/FLAG"
  setwd(master.dir)

# Load required packages
  library(ggplot2)
  library(Hmisc)
	library(doBy)
	library(plyr)

# Source functions
	source("E:/Hughes/functions_utility.r")

# Load NCA functions
  Cmax <- function(x, y) {
  # computes the Cmax
    cmax <- max(y, na.rm = T)
    tindex <- which(y == cmax)
    tmax <- x[tindex]
    unique(cmax)[1] #as there can be 2 or more equal Cmax's, choose the first
  }

  tmax <- function(x, y) {
  # computes the time of Cmax
    cmax <- max(y,na.rm=T)
    tindex <- which(y==cmax)
    tmax <- x[tindex]
    head(tmax, n=1)   #as there can be 2 or more equal Cmax's, choose the first
  }

  AUCtrapz2 <- function(x, y) {
    # computes the integral of y with respect to x using trapezoidal integration
    # does not handle missing data
    n <- length(x)
    idx <- 2:n
    AUC0t <- (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
    data.frame(AUC0t, n)
  }

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
	runname <- "RUN015_DES_1C7TAP_IBW_PPV_CORCLV_KTR_VPC"
  # processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data <- SIM.data[SIM.data$MDV == 0, ]
  if (runname == "RUN028_CELGENE") {
    SIM.data$DV <- exp(SIM.data$DV)
    SIM.data$PRED <- exp(SIM.data$PRED)
    SIM.data$IPRED <- exp(SIM.data$IPRED)
  }

# Read in the original data
  ORG.data <- read.csv("nmprep_flagged.csv", stringsAsFactors = F, na.strings = ".")
	names(ORG.data)[names(ORG.data) == "X.ID"] <- "ID"
  ORG.data <- ORG.data[ORG.data$MDV == 0 & ORG.data$FLAG == 0, ]

# -----------------------------------------------------------------------------
# Fix OCC data for STUDY == 8056
  ORG.data$DXCAT <- ORG.data$DXCATNUM
  ORG.data$OCC <- 1
  ORG.data$OCC[ORG.data$TIME >= 24 & ORG.data$TAD != 24] <- 2
  ORG.data$OCC[ORG.data$TIME >= 48 & ORG.data$TAD == 24] <- 2
  ORG.data$OCC[ORG.data$TIME >= 48 & ORG.data$TAD != 24] <- 3
  ORG.data$OCC[ORG.data$TIME >= 72 & ORG.data$TAD == 24] <- 3
  ORG.data$OCC[ORG.data$TIME >= 72 & ORG.data$TAD != 24] <- 4
  ORG.data$OCC[ORG.data$TIME >= 96 & ORG.data$TAD == 24] <- 4
  ORG.data$OCC[ORG.data$TIME >= 96 & ORG.data$TAD != 24] <- 5
  ORG.data$OCC[ORG.data$TIME >= 120 & ORG.data$TAD == 24] <- 5
  ORG.data$OCC[ORG.data$TIME >= 120 & ORG.data$TAD != 24] <- 6

# Add OCC data to SIM.data
  SIM.data$OCC <- rep(ORG.data$OCC, 1000)

  ORG.data <- rbind(
    ORG.data[(ORG.data$STUDY == 6003 | ORG.data$STUDY == 5115) & ORG.data$OCC == 1,],
    ORG.data[ORG.data$STUDY == 8056 & ORG.data$OCC %in% c(2, 3),],
    ORG.data[ORG.data$STUDY == 10016 & ORG.data$OCC %in% c(1, 5),]
  )
  SIM.data <- rbind(
    SIM.data[(SIM.data$STUDY == 6003 | SIM.data$STUDY == 5115) & SIM.data$OCC == 1,],
    SIM.data[SIM.data$STUDY == 8056 & SIM.data$OCC %in% c(2, 3),],
    SIM.data[SIM.data$STUDY == 10016 & SIM.data$OCC %in% c(1, 5),]
  )

# -----------------------------------------------------------------------------
# Calculate simulated AUC etc.
  AUCdata <- NULL
  Cmaxdata <- NULL
  NCAdata <- NULL

  AUCdata <- ddply(ORG.data, .(STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) AUCtrapz2(x$TAD, x$DV))
  Cmaxdata <- ddply(ORG.data, .(STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) Cmax(x$TAD, x$DV))
  Tmaxdata <- ddply(ORG.data, .(STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) tmax(x$TAD, x$DV))

  NCAdata.ORG <- cbind(AUCdata, "CMAX" = Cmaxdata$V1, "TMAX" = Tmaxdata$V1)

  AUCdata <- NULL
  Cmaxdata <- NULL
  NCAdata <- NULL

  AUCdata <- ddply(SIM.data, .(SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) AUCtrapz2(x$TAD, x$DV))
  Cmaxdata <- ddply(SIM.data, .(SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) Cmax(x$TAD, x$DV))
  Tmaxdata <- ddply(SIM.data, .(SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) tmax(x$TAD, x$DV))

  NCAdata.SIM <- cbind(AUCdata, "CMAX" = Cmaxdata$V1, "TMAX" = Tmaxdata$V1)

# Trim extreme values for AUC and Cmax
  trim_extreme <- function(x, lower = 0.025, upper = 0.975) {
  # Box plots mark outliers that are 1 interquartile range below and above the
  # 25th and 75th quantile
  # This function sets values outside the lower and upper quantiles to NA to
  # remove extreme outliers from boxplots
  # lower <- q25 - (q75 - q25)
  # upper <- q75 + (q75 - q25)
    lowerlimit <- quantile(x, probs = lower, na.rm = T, names = F)
    upperlimit <- quantile(x, probs = upper, na.rm = T, names = F)
    x[x < lowerlimit] <- NA
    x[x > upperlimit] <- NA
    x
  }

  NCAdata.SIM$AUC0t <- trim_extreme(NCAdata.SIM$AUC0t)
  NCAdata.SIM$CMAX <- trim_extreme(NCAdata.SIM$CMAX)

  NCAdata.SIM2 <- cbind(SOURCE = "Simulated", NCAdata.SIM)
  names(NCAdata.SIM2)

  NCAdata.ORG2 <- cbind(SOURCE = "Observed", SIM = 0, NCAdata.ORG)
  names(NCAdata.ORG2)

  NCAdata.ALL <- rbind(NCAdata.SIM2, NCAdata.ORG2)

# -----------------------------------------------------------------------------
# Assign factors to covariates
  #Assign some factors
	NCAdata.ALL$SEXf <- factor(NCAdata.ALL$SEX, labels = c("F", "M"))
	NCAdata.ALL$DXCATf <- factor(NCAdata.ALL$DXCAT, labels = c("CLL", "AML", "ALL", "MM"))
	NCAdata.ALL$STUDYf <- factor(NCAdata.ALL$STUDY)
	NCAdata.ALL$CRCLf <- factor(ifelse(NCAdata.ALL$CRCL2 <= 60, 1, 2),
	  labels = c("CrCl <60mL/min", "CrCl >60mL/min"))

  logscaleticks2 <- sort(c(
    seq(from = -10, to = 10, by = 1),
    seq(from = log10(0.000003), to = log10(300000), by = 1)
  ))

# Plot AUC
  p <- NULL
  p <- ggplot(data = NCAdata.ALL)
  p <- p + geom_boxplot(aes(x = SOURCE, y = AUC0t), notch = F)
  # p <- p + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p <- p + scale_x_discrete("Data Source")
  p <- p + scale_y_log10("AUC (ng/ml * h)")  #, lim = c(0,35)
  p <- p + ggtitle("VPC of lenalidomide exposure by AUC0-24\n")
  p <- p + facet_wrap(~STUDYf, ncol = 2, scales = "free_y")
  # p
  ggsave(paste0(runname, ".nm7/auc_ncavpc_STUDY.png"),
    width = 32, height = 42, units = c("cm"))

# Plot Cmax
  p <- NULL
  p <- ggplot(data = NCAdata.ALL)
  p <- p + geom_boxplot(aes(x = SOURCE, y = CMAX), notch = F)
  # p <- p + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p <- p + scale_x_discrete("Data Source")
  p <- p + scale_y_log10("AUC (ng/ml * h)")  #, lim = c(0,35)
  p <- p + ggtitle("VPC of lenalidomide Cmax\n")
  p <- p + facet_wrap(~STUDYf, ncol = 2, scales = "free_y")
  # p
  ggsave(paste0(runname, ".nm7/cmax_ncavpc_STUDY.png"),
    width = 32, height = 42, units = c("cm"))

# Plot Tmax
  p <- NULL
  p <- ggplot(data = NCAdata.ALL)
  p <- p + geom_boxplot(aes(x = SOURCE, y = TMAX), notch = F)
  # p <- p + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p <- p + scale_x_discrete("Data Source")
  p <- p + scale_y_log10("AUC (ng/ml * h)")  #, lim = c(0,35)
  p <- p + ggtitle("VPC of lenalidomide Cmax\n")
  p <- p + facet_wrap(~STUDYf, ncol = 2, scales = "free_y")
  # p
  ggsave(paste0(runname, ".nm7/tmax_ncavpc_STUDY.png"),
    width = 32, height = 42, units = c("cm"))
