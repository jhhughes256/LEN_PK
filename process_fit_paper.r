# Diagnostic plots for ESM material for Lena PopPK Paper
# Read data for final run: COV15/RUN016

# Remove all current objects in the workspace
	rm(list = ls(all = TRUE))
	graphics.off()

# Set the working directory
   master.dir <- "E:/Hughes/Data/PK/FLAG"
   scriptname <- "process_fit"
   setwd(master.dir)

 #Load libraries
  library(ggplot2)
	library(cowplot)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

#Customize ggplot2 theme
	theme_bw2 <- theme_set(theme_bw())
	theme_bw2 <- theme_update(
		plot.margin = unit(c(0.5, 0.5, 1, 0.5), "lines")
	)

# Read in the fit file
# Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in,".csv", sep = "")

	runfolder <- base.path # picks up the folder of the run being analysed

#Read *.fit file and attach, so column names are available
  fitdata <- read.table(
		file = file.name.in, sep = "", skip = 1, header = T,
		na.strings = c("NA","***********","1.#INFE+00")
	)

# Write to file
  write.csv(fitdata, file = file.name.out)

# Read in observed data
	nmprep <- read.csv("E:/Hughes/Data/PK/FLAG/nmprep_flagged.csv")

# Perform data specific corrections
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove dose events & missing values
  fitdata <- subset(fitdata, MDV == 0)
	fitdata$HT[fitdata$HT == 0 & fitdata$GEND == 1] <- 1.75
	fitdata$HT[fitdata$HT == 0 & fitdata$GEND == 0] <- 1.6

#Set factors and levels
	fitdata$DXCATf <- factor(fitdata$DXCAT)
	levels(fitdata$DXCATf) <- c("CLL", "AML", "ALL", "MM")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Diagnostic plots using cowplot
# Define colourblind palette and custom palette
  cbPalette <- data.frame(
		grey = as.character("#999999"),
		orange = "#E69F00",
		skyblue = "#56B4E9",
		green = "#009E73",
		yellow = "#F0E442",
		blue = "#0072B2",
		red = "#D55E00",
		pink = "#CC79A7",
		stringsAsFactors = F
	)
  myPalette <- with(cbPalette, c(blue, green, red, pink))

# Plot 1 - OBS vs PRED
# Note: legend is extracted here
  xy.lim1 <- c(-0.1, max(c(fitdata$PRED, fitdata$DV), na.rm = T))

  p1 <- NULL
  p1 <- ggplot(fitdata)
  p1 <- p1 + geom_point(aes(x = PRED, y = DV, colour = DXCATf),
	  shape = 1, alpha = 0.5, size = 1)
  p1 <- p1 + geom_abline(aes(x = PRED, y = DV),
	  intercept = 0, slope = 1, colour = "black")  # add line of identity
  p1 <- p1 + geom_smooth(aes(x = PRED, y = DV),
	  method = loess, se = T, colour = "red")  # add loess smoothing line
	p1 <- p1 + xlab("Population Predicted (ug/mL)")
	p1 <- p1 + ylab("Observed (ug/mL)")
	p1 <- p1 + coord_cartesian(xlim = xy.lim1, ylim = xy.lim1) # set limits
  p1 <- p1 + scale_colour_manual(name = "Dose Level", values = myPalette)
	pLegend <- get_legend(p1)  ### extract legend before removal ###
  p1 <- p1 + theme(legend.position = "none")  # remove legend

# Plot 2 - OBS vs IPRED
  xy.lim2 <- c(-0.1, max(c(fitdata$IPRED, fitdata$DV), na.rm = T))

  p2 <- NULL
  p2 <- ggplot(fitdata)
  p2 <- p2 + geom_point(aes(x = IPRED, y = DV, colour = DXCATf),
	  shape = 1, alpha = 0.5, size = 1)
  p2 <- p2 + geom_abline(aes(x = IPRED, y = DV),
	  intercept = 0, slope = 1, colour = "black")  # add line of identity
  p2 <- p2 + geom_smooth(aes(x = IPRED, y = DV),
	  method = loess, se = T, colour = "red")  # add loess smoothing line
  p2 <- p2 + xlab("Individual Predicted (ug/mL)")
  p2 <- p2 + ylab("Observed (ug/mL)")
	p2 <- p2 + coord_cartesian(xlim = xy.lim2, ylim = xy.lim2) # set limits
  p2 <- p2 + scale_colour_manual(name = "Dose Level", values = myPalette)
  p2 <- p2 + theme(legend.position = "none")  # remove legend

# Plot 3  CWRES vs TAD
  y.lim <- max(abs(fitdata$CWRES), na.rm = T)
  if (y.lim < 0.1) y.lim <- 0.1

  p3 <- NULL
  p3 <- ggplot(fitdata)
  p3 <- p3 + geom_point(aes(x = TAD, y = CWRES, colour = DXCATf),
	  shape = 1, alpha = 0.5, size = 1)
  p3 <- p3 + geom_abline(aes(x = TAD, y = CWRES),
	  intercept = 0, slope = 0, colour = "black")  # add zero line
  p3 <- p3 + geom_smooth(aes(x = TAD, y = CWRES),
	  method = loess, se = T, colour = "red")  # add loess smoothing line
  p3 <- p3 + xlab("Time after dose (h)")
  p3 <- p3 + ylab("CWRES")
	p3 <- p3 + coord_cartesian(ylim = c(-y.lim, y.lim))
  p3 <- p3 + scale_colour_manual(name = "Dose Level", values = myPalette)
  p3 <- p3 + theme(legend.position = "none")  # remove legend

# Plot 4   CWRES vs PRED
  p4 <- NULL
  p4 <-  ggplot(fitdata)
  p4 <- p4 + geom_point(aes(x = PRED, y = CWRES, colour = DXCATf),
	  shape = 1, alpha = 0.5, size = 1)
  p4 <- p4 + geom_abline(aes(x = PRED, y = CWRES),
	  intercept = 0, slope = 0, colour = "black")  #Add zero line
  p4 <- p4 + geom_smooth(aes(x = PRED, y = CWRES),
	  method = loess, se = T, colour = "red")        #Add loess smoothing line
  p4 <- p4 + xlab("Population Predicted (ug/mL)")
  p4 <- p4 + ylab("CWRES")
	p4 <- p4 + coord_cartesian(ylim = c(-y.lim, y.lim))
  p4 <- p4 + scale_colour_manual(name = "Dose Level", values = myPalette)
  p4 <- p4 + theme(legend.position = "none")  # remove legend

# Combine the plots with cowplot
# Define top half: two columns, align size of plots
  pTop <- plot_grid(p1, p2, ncol = 2, align = "hv")

# Define bottom half: 2 rows plots, legend on side, align size of plots
  pCWRES <- plot_grid(p3, p4, ncol = 1, align = "hv")
	pBot <- plot_grid(pCWRES, pLegend, ncol = 2, rel_widths = c(0.8, 0.2))

# Final plot: 1 column, do not alignment of plots
  plot_grid(pTop, pBot, ncol = 1)

# Save plot to file
# Cowplot implementation of ggsave, use ggplot2::ggsave() if misbehaving
	ggsave(paste(base.path, "diag_dash.png", sep = "/"),
		width = 13.05, height = 17.4, units = c("cm"))
	ggsave(paste(base.path, "diag_dash.eps", sep = "/"),
		dpi = 1200, device = cairo_ps, fallback_resolution = 1200,
	  width = 13.05, height = 17.4, units = c("cm"))
