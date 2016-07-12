###sim_check.r
##Goal: To check simulations to determine initial estimates
##NOTE: Save each plot as a diff name

# Remove any previous objects in the workspace
   rm(list=ls(all=TRUE))
   graphics.off()
   
# Set the working directory
   master.dir <- "D:/Hughes/Data"
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

# Source utility functions file
   source("D:/Hughes/functions_utility.r")
   
# Customize ggplot2 theme - R 2.15.3
   setthemebw2.1()
   
# Organise working and output directories
   working.dir <- paste(master.dir,"RAW_Clinical",sep="/")
   workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

   output.dir <- paste(working.dir,"/",scriptname,"_Output",sep="")
   if(!file.exists(output.dir)){
		dir.create(output.dir)
    }
	
   fit.file <- "D:/Hughes/Data/RAW_Clinical/nmprep_clin_Output/06003LEN_1comp1abs_SIM.nm7/06003len_1comp1abs_sim.fit"
   fitdata <- read.table(file=fit.file, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
   fitdata$Y[fitdata$Y==0] <- NA
   fitdata$DOSELVLf <- as.factor(fitdata$DOSELVL)
   
      #Conc vs TIME Week 1 per ID
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  filename.out <- paste(output.dir,"SIM1_CL12_V80_KA5",sep="/")
  to.png(plotobj,filename.out) 
  
     #Conc vs TIME Week 1 per ID (First Half)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID<36))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  to.png(plotobj,paste(filename.out,"p1"))
  
     #Conc vs TIME Week 1 per ID (Second Half)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID>35))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=DV, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  plotobj

  to.png(plotobj,paste(filename.out,"p2"))
  
  plotdata <- read.csv(paste(working.dir,"datacheck_clin_Output","plotdata.csv",sep="/"))
  
  