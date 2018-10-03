# Prepare dataset for simulation
# -----------------------------------------------------------------------------
# Dataset will then be simulated using my model and the Guglieri-Lopez model
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
  library(plyr)
  library(reshape2)  #melt

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
  theme_bw2 <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"),
    axis.title.x=element_text(size = 18, vjust = 0),
    axis.title.y=element_text(size = 18, vjust = 0, angle = 90),
    strip.text.x=element_text(size = 16),
    strip.text.y=element_text(size = 16, angle = 90),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12))

# Confidence intervals - from function utility
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)

  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)

# -----------------------------------------------------------------------------
# Read in data for plotting
# Process the simulated *.fit files
  setwd("E:/Hughes/Data/PK/FLAG/COV15")
	runname1 <- "RUN016_CL_CRCL2_FFM_VPC"
  # processSIMdata(paste(runname1,".ctl",sep=""))
  SIM.data1 <- read.csv(paste(runname1, ".nm7/", runname1, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data1 <- SIM.data1[SIM.data1$MDV == 0, ]
  SIM.data1$MODEL <- "Present Model"

  setwd(master.dir)
  runname2 <- "RUN028_CELGENE"
  # processSIMdata(paste(runname2,".ctl",sep=""))
  SIM.data2 <- read.csv(paste(runname2, ".nm7/", runname2, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data2 <- SIM.data2[SIM.data2$MDV == 0, ]
  SIM.data2$DV <- exp(SIM.data2$DV)
  SIM.data2$MODEL <- "Guglieri-Lopez et al. (2017)"

  runname3 <- "RUN029_LOPEZ"
  # processSIMdata(paste(runname3,".ctl",sep=""))
  SIM.data3 <- read.csv(paste(runname3, ".nm7/", runname3, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data3 <- SIM.data3[SIM.data3$MDV == 0, ]
  SIM.data3$MODEL <- "Connarn et al. (2017)"


# Read in the original data
  ORG.data <- read.csv("nmprep_flagged.csv", stringsAsFactors = F, na.strings = ".")
	names(ORG.data)[names(ORG.data) == "X.ID"] <- "ID"
  ORG.data <- ORG.data[ORG.data$MDV == 0 & ORG.data$FLAG == 0, ]

# Bind the simulated data by its rows
  SIM.data <- rbind(SIM.data1, SIM.data2, SIM.data3)

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
  SIM.data$OCC <- rep(ORG.data$OCC, 3000)

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

  AUCdata <- ddply(SIM.data, .(MODEL, SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) AUCtrapz2(x$TAD, x$DV))
  Cmaxdata <- ddply(SIM.data, .(MODEL, SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) Cmax(x$TAD, x$DV))
  Tmaxdata <- ddply(SIM.data, .(MODEL, SIM, STUDY, ID, OCC, SEX, CRCL2, DXCAT), function(x) tmax(x$TAD, x$DV))

  NCAdata.SIM <- cbind(AUCdata, "CMAX" = Cmaxdata$V1, "TMAX" = Tmaxdata$V1)

# Trim extreme values for AUC and Cmax
  # trim_extreme <- function(x, lower = 0.025, upper = 0.975) {
  # # Box plots mark outliers that are 1 interquartile range below and above the
  # # 25th and 75th quantile
  # # This function sets values outside the lower and upper quantiles to NA to
  # # remove extreme outliers from boxplots
  # # lower <- q25 - (q75 - q25)
  # # upper <- q75 + (q75 - q25)
  #   lowerlimit <- quantile(x, probs = lower, na.rm = T, names = F)
  #   upperlimit <- quantile(x, probs = upper, na.rm = T, names = F)
  #   x[x < lowerlimit] <- NA
  #   x[x > upperlimit] <- NA
  #   x
  # }
  #
  # NCAdata.SIM$AUC0t <- trim_extreme(NCAdata.SIM$AUC0t)
  # NCAdata.SIM$CMAX <- trim_extreme(NCAdata.SIM$CMAX)

  NCAdata.SIM2 <- cbind(SOURCE = "Sim", NCAdata.SIM)
  names(NCAdata.SIM2)

  NCAdata.ORG2 <- cbind(SOURCE = "Obs", MODEL = "None", SIM = 0, NCAdata.ORG)
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

  NCAdata.MELT <- melt(NCAdata.ALL, measure.vars = c("AUC0t", "CMAX", "TMAX"),
     value.name = "VALUE", variable.name = "METRIC")

  NCAdata.MELT$METRICf <- factor(NCAdata.MELT$METRIC)
  levels(NCAdata.MELT$METRICf) <- c("AUC 0-tlast (ng/ml * h)", "Cmax (ng/ml)", "tmax (h)")

  sumstat.ALL <- ddply(NCAdata.MELT, .(MODEL, METRICf), function(x) {
    data.frame(
      y05 = quantile(x$VALUE, 0.05, na.rm = T),
      y25 = quantile(x$VALUE, 0.25, na.rm = T),
      y50 = median(x$VALUE, na.rm = T),
      y75 = quantile(x$VALUE, 0.75, na.rm = T),
      y95 = quantile(x$VALUE, 0.95, na.rm = T)
    )
  })

# Define colourblind palette and custom palette
  cbPalette <- data.frame(
		grey = "#999999",
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  myPalette <- with(cbPalette, c(green, red, blue, pink))

# Create plot comparing NCA metrics between models
  p <- NULL
  p <- ggplot()
  p <- p + geom_point(aes(x = MODEL, y = y50, colour = MODEL),
    data = sumstat.ALL[sumstat.ALL$MODEL != "None", ], size = 2.8)
  p <- p + geom_errorbar(aes(x = MODEL, ymin = y05, ymax = y95, colour = MODEL),
    data = sumstat.ALL[sumstat.ALL$MODEL != "None", ], size = 1.2)
  p <- p + geom_hline(aes(yintercept = y50),
    data = sumstat.ALL[sumstat.ALL$MODEL == "None", ], size = 1)
  p <- p + geom_hline(aes(yintercept = y05), linetype = "dashed",
    data = sumstat.ALL[sumstat.ALL$MODEL == "None", ], size = 1)
  p <- p + geom_hline(aes(yintercept = y95), linetype = "dashed",
    data = sumstat.ALL[sumstat.ALL$MODEL == "None", ], size = 1)
  p <- p + facet_wrap(~METRICf, nrow = 3, scales = "free_x")
  p <- p + scale_colour_manual(values = myPalette)
  p <- p + guides(colour = F)
  p <- p + xlab(NULL)
  p <- p + ylab(NULL)
  p <- p + coord_flip()
  p

# Save to file
  ggsave("ncavpc_compare_v1.png",
    width = 10, height = 12.78, units = c("in"), dpi = 72)
  ggsave("ncavpc_compare_v1.eps",
    dpi = 1200, device = cairo_ps, fallback_resolution = 1200,
    width = 17.4, height = 23.4, units = c("cm"))
