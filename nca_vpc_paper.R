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

# Plot AUC
  p1 <- NULL
  p1 <- ggplot(data = NCAdata.ALL)
  p1 <- p1 + geom_boxplot(aes(x = SOURCE, y = AUC0t), notch = F)
  # p1 <- p1 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p1 <- p1 + scale_x_discrete("\n")
  p1 <- p1 + scale_y_log10("AUC (ng/ml * h)\n", breaks = logscaleticks2)  #, lim = c(0,35)

# Plot Cmax
  p2 <- NULL
  p2 <- ggplot(data = NCAdata.ALL)
  p2 <- p2 + geom_boxplot(aes(x = SOURCE, y = CMAX), notch = F)
  # p2 <- p2 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p2 <- p2 + scale_x_discrete("\nData Source")
  p2 <- p2 + scale_y_log10("Cmax (ng/ml)\n", breaks = logscaleticks2)  #, lim = c(0,35)

# Plot Tmax
  p3 <- NULL
  p3 <- ggplot(data = NCAdata.ALL)
  p3 <- p3 + geom_boxplot(aes(x = SOURCE, y = TMAX), notch = F)
  # p3 <- p3 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p3 <- p3 + scale_x_discrete("\n")
  p3 <- p3 + scale_y_log10("tmax (ng/ml * h)\n", breaks = logscaleticks2)  #, lim = c(0,35)

  p4 <- plot_grid(p1, p2, p3, nrow = 1)
  p4
  ggsave(paste0(runname, ".nm7/ncavpc_ALL.png"),
    width = 32, height = 21, units = c("cm"))

# Identity boxplots
  p1 <- NULL
  p1 <- ggplot(data = sumstat.ALL)
  p1 <- p1 + geom_boxplot(aes(x = SOURCE, ymin = AUCy05, lower = AUCy25, middle = AUCy50, upper = AUCy75, ymax = AUCy95),
    notch = F, stat = "identity")
  # p1 <- p1 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p1 <- p1 + scale_x_discrete("\n")
  p1 <- p1 + scale_y_log10("AUC (ng/ml * h)", breaks = logscaleticks2)  #, lim = c(0,35)

  # Plot Cmax
  p2 <- NULL
  p2 <- ggplot(data = sumstat.ALL)
  p2 <- p2 + geom_boxplot(aes(x = SOURCE, ymin = CMAXy05, lower = CMAXy25, middle = CMAXy50, upper = CMAXy75, ymax = CMAXy95),
    notch = F, stat = "identity")
  # p2 <- p2 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p2 <- p2 + scale_x_discrete("\nData Source          ")
  p2 <- p2 + scale_y_log10("Cmax (ng/ml)", breaks = logscaleticks2[-0.03])  #, lim = c(0,35)

  # Plot Tmax
  p3 <- NULL
  p3 <- ggplot(data = sumstat.ALL)
  p3 <- p3 + geom_boxplot(aes(x = SOURCE, ymin = TMAXy05, lower = TMAXy25, middle = TMAXy50, upper = TMAXy75, ymax = TMAXy95),
    notch = F, stat = "identity")
  # p3 <- p3 + stat_summary(aes_string(x = covname, y = "CL"), data = statdata,
  #    fun.data = boxplot.give.n, geom = "text", size = 6, colour = "red")
  p3 <- p3 + scale_x_discrete("\n")
  p3 <- p3 + scale_y_log10("tmax (hours)", breaks = logscaleticks2)  #, lim = c(0,35)

  p4 <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1.01, 1.08, 0.899),
    labels = c("A", "B", "C"), label_fontface = "plain")
  p4
  ggsave(paste0(runname, ".nm7/ncavpc_ALL_identity.png"),
    width = 17.4, height = 15.7, units = c("cm"))
  ggsave(paste0(runname, ".nm7/ncavpc_ALL_identity.eps"),
    dpi = 1200, device = cairo_ps, fallback_resolution = 1200,
    width = 17.4, height = 15.7, units = c("cm"))


# Plot Predicted vs. Observed Metrics
  NCAdata.SIM3 <- melt(NCAdata.SIM2, measure.vars = c("AUC0t", "CMAX", "TMAX"), value.name = "METRIC.SIM")
  NCAdata.OBS3 <- melt(NCAdata.ORG2, measure.vars = c("AUC0t", "CMAX", "TMAX"), value.name = "METRIC.OBS")

  NCAdata.MELT2 <- ddply(NCAdata.SIM3, .(SIM), function(x) {
    x$METRIC.OBS <- NCAdata.OBS3$METRIC.OBS
    x
  })

  p5 <- NULL
  p5 <- ggplot(data = NCAdata.MELT2)
  p5 <- p5 + geom_point(aes(x = METRIC.OBS, y = METRIC.SIM), alpha = 0.5)
  p5
