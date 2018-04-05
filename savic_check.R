### Check if transits are working via NONMEM simulation
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data/PK/CHECK_SAVIC"
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
	runname <- "RUN001_TRANSIT_SAVIC"
  # processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data1 <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data1 <- SIM.data1[SIM.data1$EVID == 0, ]
  SIM.data1$SET <- "SAVIC"

  runname <- "RUN002_TRANSIT_MULTI"
  processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data2 <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data2 <- SIM.data2[SIM.data2$EVID == 0, ]
  SIM.data2$SET <- "SUPER"

  runname <- "RUN003_HARDCODED"
  # processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data3 <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data3 <- SIM.data3[SIM.data3$EVID == 0, ]
  SIM.data3$SET <- "HARDCODED"

  runname <- "RUN004_TRANSIT_SAVIC_WILKINS"
  processSIMdata(paste(runname,".ctl",sep=""))
  SIM.data4 <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  SIM.data4 <- SIM.data4[SIM.data4$EVID == 0, ]
  SIM.data4$SET <- "WILKINS"

  SIM.data <- rbind(SIM.data1, SIM.data2, SIM.data3, SIM.data4)

# -----------------------------------------------------------------------------
# Plot data
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = TIME, y = CD, group = ID), data = SIM.data[SIM.data$ID == 3, ],
	  geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.3)
	p <- p + stat_summary(aes(x = TIME, y = CD, group = ID), data = SIM.data[SIM.data$ID == 3, ],
		fun.y = median, geom = "line", colour = "black", size = 1)
  p + facet_wrap(~ID+SET, ncol = 4)

  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = TIME, y = DV, group = ID), data = SIM.data[SIM.data$ID == 3, ],
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.3)
  p <- p + stat_summary(aes(x = TIME, y = DV, group = ID), data = SIM.data[SIM.data$ID == 3, ],
    fun.y = median, geom = "line", colour = "black", size = 1)
  p + facet_wrap(~ID+SET, ncol = 4)

# -----------------------------------------------------------------------------
# Check original code
  runname <- "w14_sim"
  processSIMdata(paste(runname,".ctl",sep=""))
  ORG.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  ORG.data <- ORG.data[ORG.data$EVID == 0, ]

  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = TIME, y = CD), data = ORG.data,
		fun.y = median, geom = "line", colour = "black", size = 1)
  p

# -----------------------------------------------------------------------------
# Check TAD used in Savic Wilkins method
  runname <- "RUN001_TRANSIT_SAVIC_TAD"
  processSIMdata(paste(runname,".ctl",sep=""))
  TAD.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  TAD.data <- TAD.data[TAD.data$EVID == 0, ]

  p <- NULL
  p <- ggplot(data = TAD.data[TAD.data$ID == 3, ])
  p <- p + geom_line(aes(x=TIME, y=TAD))
  p

# -----------------------------------------------------------------------------
# Check NDOSE and I in Superimposition methods
  runname <- "RUN002_TRANSIT_MULTI_NDOSE"
  processSIMdata(paste(runname,".ctl",sep=""))
  NDI.data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  NDI.data <- NDI.data[NDI.data$EVID == 0, ]

  p <- NULL
  p <- ggplot(data = NDI.data)
  p <- p + stat_summary(aes(x = TIME, y = CD), data = NDI.data,
		fun.y = median, geom = "line", colour = "black", size = 1)
  p + facet_wrap(~ID, ncol = 3)
