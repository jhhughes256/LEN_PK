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
  library(plyr)

# Source functions
	source("E:/Hughes/functions_utility.r")

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
# Load the three separate simulation files
	# runname <- "RUN003_BASE15"
  # runname <- "RUN006_BASE50"
  runname <- "RUN009_BASE15000"
  # processSIMdata(paste(runname,".ctl",sep=""))
  base_data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  base_data <- base_data[base_data$MDV == 0, ]

  # runname <- "RUN004_LOPEZ15"
  # runname <- "RUN007_LOPEZ50"
  runname <- "RUN010_LOPEZ15000"
  # processSIMdata(paste(runname,".ctl",sep=""))
  lopez_data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  lopez_data <- lopez_data[lopez_data$MDV == 0, ]

  # runname <- "RUN005_CELGENE15"
  # runname <- "RUN008_CELGENE50"
  runname <- "RUN011_CELGENE15000"
  # processSIMdata(paste(runname,".ctl",sep=""))
  celgene_data <- read.csv(paste(runname, ".nm7/", runname, ".fit.csv", sep = ""),
    stringsAsFactors = F, na.strings = ".")
  celgene_data <- celgene_data[celgene_data$MDV == 0, ]
  celgene_data$DV <- exp(celgene_data$DV)

# -----------------------------------------------------------------------------
# Time binning
  timebin_cuts <- c(0.52, 1.02, 1.52, 2.02, 2.52, 3.02, 3.52, 4.02, 5.02, 6.02,
    7.02, 8.02, 10.02, 12.02, 14.02, 16.02, 18.02, 20.02, 22.02, 24.02)
  base_data$TIMEBIN <- cut2(base_data$TIME, cuts = timebin_cuts, levels.mean = T)
  base_data$TIMEBIN <- as.numeric(paste(base_data$TIMEBIN))

  lopez_data$TIMEBIN <- cut2(lopez_data$TIME, cuts = timebin_cuts, levels.mean = T)
  lopez_data$TIMEBIN <- as.numeric(paste(lopez_data$TIMEBIN))

  celgene_data$TIMEBIN <- cut2(celgene_data$TIME, cuts = timebin_cuts, levels.mean = T)
  celgene_data$TIMEBIN <- as.numeric(paste(celgene_data$TIMEBIN))

# Covariates
# Assign factors to covariates
  base_data$IDf <- as.factor(base_data$ID)
  base_data$SEXf <- factor(base_data$SEX, labels=c("F", "M"))
  base_data$CRCLf <- factor(ifelse(base_data$CRCL2<=60,1,2),
    labels=c("CrCl <60mL/min","CrCl >60mL/min"))

  lopez_data$IDf <- as.factor(lopez_data$ID)
  lopez_data$SEXf <- factor(lopez_data$SEX, labels=c("F", "M"))
  lopez_data$CRCLf <- factor(ifelse(lopez_data$CRCL2<=60,1,2),
    labels=c("CrCl <60mL/min","CrCl >60mL/min"))

  celgene_data$IDf <- as.factor(celgene_data$ID)
  celgene_data$SEXf <- factor(celgene_data$SEX, labels=c("F", "M"))
  celgene_data$CRCLf <- factor(ifelse(celgene_data$CRCL2<=60,1,2),
    labels=c("CrCl <60mL/min","CrCl >60mL/min"))

# -----------------------------------------------------------------------------
# Uppsala Style VPC
# Calculate 5, 50 and 95 percentiles for each simulated study (S)
  # base_data_bysim <- ddply(base_data, .(SIM, TIMEBIN), function(x) {
  base_data_bysim <- ddply(base_data, .(STUDY, TIMEBIN), function(x) {
    data.frame(
      Model = "Present Model",
      medianS = median(x$DV*1000),
      loCI90S = CI90lo(x$DV*1000),
      hiCI90S = CI90hi(x$DV*1000)
    )
  })

  # lopez_data_bysim <- ddply(lopez_data, .(SIM, TIMEBIN), function(x) {
  lopez_data_bysim <- ddply(lopez_data, .(STUDY, TIMEBIN), function(x) {
    data.frame(
      Model = "Guglieri-Lopez et al. (2017)",
      medianS = median(x$DV*1000),
      loCI90S = CI90lo(x$DV*1000),
      hiCI90S = CI90hi(x$DV*1000)
    )
  })

  # celgene_data_bysim <- ddply(celgene_data, .(SIM, TIMEBIN), function(x) {
  celgene_data_bysim <- ddply(celgene_data, .(STUDY, TIMEBIN), function(x) {
    data.frame(
      Model = "Connarn et al. (2017)",
      medianS = median(x$DV*1000),
      loCI90S = CI90lo(x$DV*1000),
      hiCI90S = CI90hi(x$DV*1000)
    )
  })

  data_bysim <- rbind(base_data_bysim, lopez_data_bysim, celgene_data_bysim)

  p2 <- NULL

  titletext <- "VPC - Uppsala Style\n"
  p2 <- ggplot(data = data_bysim)
  # p2 <- p2 + ggtitle(titletext)

  p2 <- p2 + stat_summary(aes(x = TIMEBIN, y = medianS, fill = Model),
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.33)
  p2 <- p2 + stat_summary(aes(x = TIMEBIN, y = medianS, linetype = Model),
    fun.y = median, geom = "line", size = 1)

  p2 <- p2 + stat_summary(aes(x = TIMEBIN, y = loCI90S, linetype = Model),
    fun.y = median, geom = "line", size = 0.8, colour = "grey40")
  p2 <- p2 + stat_summary(aes(x = TIMEBIN, y = hiCI90S, linetype = Model),
    fun.y = median, geom = "line", size = 0.8, colour = "grey40")

  p2 <- p2 + scale_linetype_manual(values = c("solid", "longdash", "dotted"))
  p2 <- p2 + scale_fill_manual(values = rep("#999999", 3))
  p2 <- p2 + theme(legend.position = c(0.75, 0.85), legend.direction = "vertical")

  p2 <- p2 + scale_y_log10("Lenalidomide Concentration (mg/L)\n",
    breaks = c(1, 10, 100, 1000))
  p2 <- p2 + scale_x_continuous("\nTime (hours)", breaks = 0:6*4)
  p2 <- p2 + coord_cartesian(ylim = c(1, 1000))
  p2

  #normal scale- non-facetted
  ggsave("Uppsala_VPC_v2BW.png", width=17.4, height=14, units=c("cm"))
  ggsave("Uppsala_VPC_v2BW.eps", width=17.4, height=14, units=c("cm"),
    dpi = 1200, device = cairo_ps, fallback_resolution = 1200)

  # p <- NULL
  #
  # titletext <- "VPC - Uppsala Style\n"
  # p <- ggplot(aes(colour = Model), data = data_bysim)
  # p <- p + ggtitle(titletext)
  #
  # # p <- p + stat_summary(aes(x = TIMEBIN, y = medianS),
  # #   geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.3)
  # p <- p + stat_summary(aes(x = TIMEBIN, y = medianS),
  #   fun.y = median, geom = "line", size = 1)
  # # p <- p + stat_summary(aes(x = TIMEBIN, y = loCI90S),
  # #   geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.1)
  # p <- p + stat_summary(aes(x = TIMEBIN, y = loCI90S),
  #   fun.y = median, geom = "line", linetype = "dashed", size = 0.8)
  # # p <- p + stat_summary(aes(x = TIMEBIN, y = hiCI90S),
  # # 	geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", fill = "red", alpha = 0.1)
  # p <- p + stat_summary(aes(x = TIMEBIN, y = hiCI90S),
  #   fun.y = median, geom = "line", linetype = "dashed", size = 0.8)
  #
  # p <- p + scale_colour_manual(values = c("red", "blue", "green4"))
  # p <- p + scale_y_log10("Lenalidomide Concentration (mg/L)\n")
  # p <- p + scale_x_continuous("\nTime (hours)", breaks = 0:12*2)
  # p <- p + coord_cartesian(ylim = c(0.0001, 1))
  # p
