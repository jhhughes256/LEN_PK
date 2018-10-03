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
  library(reshape2)
  library(cowplot)
  library(scales)

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
	runname <- "RUN016_CL_CRCL2_FFM_VPC"
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

  NCAdata.SIM2 <- cbind(SOURCE = "Sim", NCAdata.SIM)
  names(NCAdata.SIM2)

  NCAdata.ORG2 <- cbind(SOURCE = "Obs", SIM = 0, NCAdata.ORG)
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

  logscaleticks2 <- 10^sort(c(
    seq(from = -10, to = 10, by = 1),
    seq(from = log10(0.000003), to = log10(300000), by = 1)
  ))

  NCAdata.MELT <- melt(NCAdata.ALL, measure.vars = c("AUC0t", "CMAX", "TMAX"))

  sumstat.ALL <- ddply(NCAdata.ALL, .(SOURCE), function(x) {
    data.frame(
      AUCy05 = quantile(x$AUC0t, 0.05, na.rm = T),
      AUCy25 = quantile(x$AUC0t, 0.25, na.rm = T),
      AUCy50 = median(x$AUC0t, na.rm = T),
      AUCy75 = quantile(x$AUC0t, 0.75, na.rm = T),
      AUCy95 = quantile(x$AUC0t, 0.95, na.rm = T),
      CMAXy05 = quantile(x$CMAX, 0.05, na.rm = T),
      CMAXy25 = quantile(x$CMAX, 0.25, na.rm = T),
      CMAXy50 = median(x$CMAX, na.rm = T),
      CMAXy75 = quantile(x$CMAX, 0.75, na.rm = T),
      CMAXy95 = quantile(x$CMAX, 0.95, na.rm = T),
      TMAXy05 = quantile(x$TMAX, 0.05, na.rm = T),
      TMAXy25 = quantile(x$TMAX, 0.25, na.rm = T),
      TMAXy50 = median(x$TMAX, na.rm = T),
      TMAXy75 = quantile(x$TMAX, 0.75, na.rm = T),
      TMAXy95 = quantile(x$TMAX, 0.95, na.rm = T)
    )
  })

# Define colourblind palette
  cbPalette <- c("#000000", "#999999")

# Plot using density plots?
  # p1 <- NULL
  # p1 <- ggplot()
  # p1 <- p1 + geom_density(aes(x = AUC0t, y = ..density.., colour = SOURCE),
  #   data = NCAdata.ALL)
  # p1 <- p1 + geom_vline(aes(xintercept = AUCy50, colour = SOURCE),
  #   data = sumstat.ALL)
  # p1 <- p1 + geom_vline(aes(xintercept = AUCy05, colour = SOURCE),
  #   data = sumstat.ALL)
  # p1 <- p1 + geom_vline(aes(xintercept = AUCy95, colour = SOURCE),
  #   data = sumstat.ALL)
  # p1 <- p1 + scale_x_log10("AUC (ng/ml * h)\n", breaks = logscaleticks2)  #, lim = c(0,35)
  # p1 <- p1 + ylab("Density")
  # p1 <- p1 + scale_colour_manual(name = "Data", values = cbPalette)

# Plot using histogram
# AUC
  p1 <- NULL
  p1 <- ggplot()
  p1 <- p1 + stat_bin(aes(x = AUC0t, y = ..count../sum(..count..)), bins = 12,
    data = NCAdata.ALL[NCAdata.ALL$SOURCE == "Sim", ], colour = "black", fill = "white")
  p1 <- p1 + geom_vline(aes(xintercept = AUCy50, colour = SOURCE),
    data = sumstat.ALL, size = 1.2)
  p1 <- p1 + geom_vline(aes(xintercept = AUCy05, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p1 <- p1 + geom_vline(aes(xintercept = AUCy95, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p1 <- p1 + scale_x_log10("AUC (ng/ml * h)", breaks = logscaleticks2)  #, lim = c(0,35)
  p1 <- p1 + scale_y_continuous("", labels = percent, expand = c(0, 0), limits = c(0, 0.181), breaks = 0:4*0.04)
  p1 <- p1 + scale_colour_manual(name = "Data", values = cbPalette)
  p1

# Cmax
  p2 <- NULL
  p2 <- ggplot()
  p2 <- p2 + stat_bin(aes(x = CMAX, y = ..count../sum(..count..)), bins = 12,
    data = NCAdata.ALL[NCAdata.ALL$SOURCE == "Sim", ], colour = "black", fill = "white")
  p2 <- p2 + geom_vline(aes(xintercept = CMAXy50, colour = SOURCE),
    data = sumstat.ALL, size = 1.2)
  p2 <- p2 + geom_vline(aes(xintercept = CMAXy05, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p2 <- p2 + geom_vline(aes(xintercept = CMAXy95, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p2 <- p2 + scale_x_log10("Cmax (ng/ml)", breaks = logscaleticks2)  #, lim = c(0,35)
  p2 <- p2 + scale_y_continuous("Percentage Count\n", labels = percent, expand = c(0, 0), limits = c(0, 0.167))

  p2 <- p2 + scale_colour_manual(name = "Data", values = cbPalette)
  p2

# tmax
  p3 <- NULL
  p3 <- ggplot()
  p3 <- p3 + stat_bin(aes(x = TMAX, y = ..count../sum(..count..)), bins = 12,
    data = NCAdata.ALL[NCAdata.ALL$SOURCE == "Sim", ], colour = "black", fill = "white")
  p3 <- p3 + geom_vline(aes(xintercept = TMAXy50, colour = SOURCE),
    data = sumstat.ALL, size = 1.2)
  p3 <- p3 + geom_vline(aes(xintercept = TMAXy05, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p3 <- p3 + geom_vline(aes(xintercept = TMAXy95, colour = SOURCE),
    data = sumstat.ALL, linetype = "dashed", size = 1.2)
  p3 <- p3 + scale_x_log10("tmax (hours)", breaks = logscaleticks2)  #, lim = c(0,35)
  p3 <- p3 + scale_y_continuous("", labels = percent, expand = c(0, 0), limits = c(0, 0.202), breaks = 0:5*0.04)
  p3 <- p3 + scale_colour_manual(name = "Data", values = cbPalette)
  p3

  p4 <- plot_grid(p1, p2, p3, ncol = 1, align = "vh", hjust = -3.2,
    labels = c("A", "B", "C"))
  p4
  ggsave(paste0(runname, ".nm7/ncavpc_ALL_identity_v2_BW.png"),
    width = 17.4, height = 23.4, units = c("cm"))
  ggsave(paste0(runname, ".nm7/ncavpc_ALL_identity_v2_BW.eps"),
    dpi = 1200, device = cairo_ps, fallback_resolution = 1200,
    width = 17.4, height = 23.4, units = c("cm"))
