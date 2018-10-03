###nmprep.r
##Goal: To collate tables of missing data contained within nonclinical raw data obtained on 23rd March 2016
##Note: Based heavily off of datacheck_cyt_script2.r -> Richards code

# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master.dir <- "E:/Hughes/Data"
  scriptname <- "nmprep_clin"
  setwd(master.dir)

# Load libraries
  library(ggplot2)
  library(doBy)
  library(Hmisc)
  library(plyr)
  library(grid)
  library(reshape)
  library(stringr)
  library(scales)
  library(cowplot)
  library(gridExtra)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

# Customize ggplot2 theme - R 2.15.3
  theme_bw2 <- theme_set(theme_bw(base_size = 22))
  theme_bw2 <- theme_update(plot.margin = unit(c(1, 0.5, 3, 0.5), "lines"),
    axis.title.x = element_text(size = 18, vjust = 0),
    axis.title.y = element_text(size = 18, vjust = 0, angle = 90),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16, angle = 90),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16))

# Organise working and output directories
  working.dir <- paste(master.dir,"RAW_Clinical",sep="/")
  workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

  output.dir <- paste(working.dir,"/",scriptname,"_Output",sep="")
  if(!file.exists(output.dir)){
	dir.create(output.dir)
  }

  filename <- paste(output.dir,"nmprep_flagged.csv",sep="/")
  nmprep <- read.csv(filename, na.strings = ".")


  locf <- function (x) {
  #Last observation carried forward
  #Finds an NA and carries forward the previous value
    good <- !is.na(x)
    positions <- seq(length(x))
    good.positions <- good * positions
    last.good.position <- cummax(good.positions)
    last.good.position[last.good.position == 0] <- NA
    x[last.good.position]
  }

  nmprep$DOSE <- locf(nmprep$AMT)
  nmprep$DVNORM <- nmprep$DV/nmprep$DOSE

  bin_cuts <- c(0.52, 1.02, 2.02, 3.02, 5.02, 8.02, 49)
  nmprep$TADBIN <- cut2(nmprep$TAD, cuts = bin_cuts, levels.mean = T)
  levels(nmprep$TADBIN)[length(bin_cuts)] <- 24
  nmprep$TADBIN <- as.numeric(paste(nmprep$TADBIN))

  dose_bins <- c(8, 26, 80)
  nmprep$DOSEf <- cut2(nmprep$DOSE, cuts = dose_bins)
  levels(nmprep$DOSEf) <- c("<10mg", "10-25mg", ">25mg")

  nmprep$DXCATf <- factor(nmprep$DXCATNUM)
  levels(nmprep$DXCATf) <- c("CLL", "AML", "ALL", "MM")

# Create plot function
  dvnormPlot <- function(xCol, guideName) {
    p <- NULL
    p <- ggplot(aes(x = TAD, y = DVNORM*1000), data = nmprep)
    p <- p + geom_point(colour = "#999999", alpha = 0.2)
    p <- p + stat_summary(aes(x = TADBIN, y = DVNORM*1000, linetype = get(xCol)), fun.y = median,
      geom = "line", size = 1.2)
    p <- p + scale_y_log10(NULL, labels = comma)
    p <- p + scale_x_continuous(NULL, breaks = 0:6*4)
    p <- p + scale_linetype_manual(guideName, values = c(1, 2, 3, 6))
    p
  }

# Create plots and use cowplot to create grid
  p1 <- dvnormPlot("DOSEf", "Dosage")
  p2 <- dvnormPlot("DXCATf", "Cancer")
  p3 <- plot_grid(p1, p2, align = "vh", labels = c("A", "B"),
    ncol = 1, hjust = -5)

# Create text grobs for common y and x axis labels
  y.grob <- textGrob("Dose Normalised Concentrations (ng/mL)\n", vjust = 0.7,
    gp = gpar(fontface = "plain", col = "black", fontsize = 18), rot = 90)
  x.grob <- textGrob("Time After Last Dose (hours)", hjust = 0.6, vjust = -1,
    gp = gpar(fontface = "plain", col = "black", fontsize = 18))

# Produce final figure
  plot_grid(grid.arrange(arrangeGrob(p3, left = y.grob, bottom = x.grob)))

  ggsave("dvnormplot_v2_BW.png", width = 17.4, height = 23.4, units = c("cm"))
  ggsave("dvnormplot_v2_BW.eps", width = 17.4, height = 23.4, units = c("cm"),
    dpi = 1200, device = cairo_ps, fallback_resolution = 1200)
