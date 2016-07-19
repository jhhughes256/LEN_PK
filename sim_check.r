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
   levels(fitdata$DOSELVLf) <- paste("Dose Level",levels(fitdata$DOSELVLf))
   
#------------------------------- IPRED vs TIME Facetted for ID (no overlay)
  filename.out <- paste(output.dir,"SIM8_CL12_V60_KA4",sep="/")

     #Conc vs TIME Week 1 per ID (1/4)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID<46))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=IPRED, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj

  to.png(plotobj,paste(filename.out,"IDfacet_p1",sep="_"))
  
     #Conc vs TIME Week 1 per ID (2/4)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID>45&ID<92))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=IPRED, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj

  to.png(plotobj,paste(filename.out,"IDfacet_p2",sep="_"))
	
     #Conc vs TIME Week 1 per ID (3/4)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID>91&ID<138))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=IPRED, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj

  to.png(plotobj,paste(filename.out,"IDfacet_p3",sep="_"))
			 
			 #Conc vs TIME Week 1 per ID (4/4)
  plotobj <- NULL
  titletext <- paste("Observed Concentrations in Week 1\n")
  plotobj <- ggplot(data=subset(fitdata,ID>137))
  plotobj <-  plotobj + geom_point(aes(x=TIME, y=IPRED, colour=DOSELVLf), size=3, alpha=0.5)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj

  to.png(plotobj,paste(filename.out,"IDfacet_p4",sep="_"))
  
  
#------------------------------- IPRED w/ OBS vs TIME Facetted for ID 	-> COLOUR = DOSELVLf
#pull observed data from datacheck_clin.r into r environment
  plotdata <- read.csv(paste(output.dir,"fulldata.csv",sep="/"), stringsAsFactors=F, na.strings=c("."))
  
  plotobj <- NULL
  titletext <- paste("log IPRED on OBSERVED - Week 1\n")
  plotobj <- ggplot(data=subset(fitdata))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=plotdata, size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
  
  to.png(plotobj,paste(filename.out,"logIPREDall",sep="_")) 
    
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 0 and 3 (25mg)\n")
  doseflag <- 0
  doseflag2 <- 3
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag2), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", data=subset(fitdata,DOSELVL==doseflag2), size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d03",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 1 (2.5 - 5mg)\n")
  doseflag <- 1
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d1",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 2 (2.5 - 7.5mg)\n")
  doseflag <- 2
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d2",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 4 35mg)\n")
  doseflag <- 4
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d4",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 5 35mg)\n")
  doseflag <- 5
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d5",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("IPRED on OBSERVED - Dose Level 6 35mg)\n")
  doseflag <- 6
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~ID)
  #plotobj
 
  to.png(plotobj,paste(filename.out,"IPRED_d6",sep="_")) 
  
#------------------------------- PRED w/ OBS vs TIME 					-> COLOUR = DOSELVLf
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - All Dose Levels\n")
  plotobj <- ggplot(data=subset(fitdata))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=plotdata, size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_log10("Concentration (ng/ml)",lim=c(0.001,5))
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  plotobj <- plotobj + facet_wrap(~DOSELVLf)
  #plotobj
  
  to.png(plotobj,paste(filename.out,"logPREDall",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 0 and 3 (25mg)\n")
  doseflag <- 0
  doseflag2 <- 3
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag2), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), data=subset(fitdata,DOSELVL==doseflag2), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d03",sep="_")) 
    
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 1 (2.5mg - 5mg)\n")
  doseflag <- 1
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d1",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 2 (2.5mg - 7.5mg)\n")
  doseflag <- 2
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d2",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 4 (35mg)\n")
  doseflag <- 4
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d4",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 5 (50mg)\n")
  doseflag <- 5
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d5",sep="_")) 
  
  plotobj <- NULL
  titletext <- paste("PRED on OBSERVED - Dose Level 6 (75mg)\n")
  doseflag <- 6
  plotobj <- ggplot(data=subset(fitdata,DOSELVL==doseflag))
  plotobj <- plotobj + geom_point(aes(x=TIME, y=DV,colour=DOSELVLf), data=subset(plotdata,DOSELVL==doseflag), size=3, alpha=0.5)
  plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour="black", size=1, alpha =1)
  plotobj <- plotobj + ggtitle(titletext) #+ theme(legend.position="none")                   
  plotobj <- plotobj +  scale_y_continuous("Concentration (ng/ml)")
  plotobj <- plotobj +  scale_x_continuous("Time after first dose (hours)", lim=c(0,24))  #, lim=c(0,60), breaks=seq(from=0, to=60, by=24)
  plotobj <- plotobj + scale_colour_discrete("Dose Level")
  #plotobj
  
  to.png(plotobj,paste(filename.out,"PRED_d6",sep="_")) 