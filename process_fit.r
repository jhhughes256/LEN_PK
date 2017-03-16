###process_fit.r
##Goal: To process fit file form nonmem for analysis of model
##Note: Based heavily off of process_fit_10July16.r -> DF code

#remove all current objects in the workspace
	rm(list=ls(all=TRUE))
	graphics.off()

# Set the working directory
   master.dir <- "E:/Hughes/Data/PK/REDO"
   scriptname <- "process_fit"
   setwd(master.dir)

 #Load libraries
  library(ggplot2)
  library(doBy)
  library(plyr)
  library(stringr)
  library(R2HTML)
  library(Hmisc)
  library(grid)
  library(reshape2)
  library(GGally)

# Source utility functions file
  source("E:/Hughes/functions_utility.r")

#Use custom ggplot 2 theme
	theme_custom <- theme_set(theme_bw(18))
	theme_custom <- theme_update(plot.margin = unit(c(1,0.5,3,0.5), "lines"))
#stops the y axis lable mashing onto the axis values			#NOTE: untested code from DF
  #theme_custom <- theme_update(axis.title.y=element_text(vjust=1, angle = 90))

#--------------------------------------------------------------------------------------------------
#Customize ggplot2 theme - R 2.15.3
	theme_bw2 <- theme_set(theme_bw(base_size = 22))
	theme_bw2 <- theme_update(plot.margin = unit(c(1,1,3,1), "lines"),
														axis.title.x=element_text(size = 16, vjust = 0),
														axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
														strip.text.x=element_text(size = 16),
														strip.text.y=element_text(size = 16, angle = 90))

#--------------------------------------------------------------------------------------------------

#Define metadata  #GOLD

#Comment
	metacomment <- NULL
	metacomment <-  "Lenalidomide PK"

#Generating script - this is added using a macro defined in Notepad++ - invoked using ctl-alt-M
	metafilepath <- NULL
	metafilepath <- "E:/Hughes/lena_scripts/process_fit.r"
	scriptname <- gsub(".r","",basename(metafilepath))
	working.dir <- dirname(metafilepath)

#Last run-time of script
	metadatetime <- NULL
	metadatetime <- Sys.time()

	metadata1 <- paste("#Comment: ",metacomment,sep="")
	metadata2 <- paste("#Generating script: ",metafilepath,sep="")
	metadata3 <- paste("#Script run: ",metadatetime,sep="")
	metadata4 <- paste("#Script run: ",metadatetime,sep="")
#--------------------------------------------------------------------------------------------------

#Read in the fit file
#Set the name of the required file and set the working directory
  cat("Select one of files in directory to process:\n")
  path <- gsub("\\\\", "/", file.choose())
  base.path <- dirname(path)
  setwd(base.path)
  file.name.in <- basename(path)
  file.name.out <- paste(file.name.in,".csv", sep="")

	runfolder <- base.path #picks up the folder of the run being analysed


#Read *.fit file and attach, so column names are available
  fitdata <- read.table(file=file.name.in, sep="", skip=1, header=T, na.strings=c("NA","***********","1.#INFE+00"))
#Write to file
  write.csv(fitdata, file=file.name.out)

	nmprep <- read.csv("E:/Hughes/Data/PK/REDO/COV24/nmprep_allstudies.csv")

#Remove dose events & missing values
  fitdata <- subset(fitdata, MDV==0)
	fitdata$HT[fitdata$HT==0&fitdata$GEND==1] <- 1.75
	fitdata$HT[fitdata$HT==0&fitdata$GEND==0] <- 1.6

#Set factors and levels (not all categorical covariates have been included!)
	fitdata$IDf <- as.factor(fitdata$ID)

	fitdata$STUDYf <- factor(fitdata$STUDY)
	levels(fitdata$STUDYf) <- paste("Study",levels(fitdata$STUDYf))

	fitdata$OCCf <- factor(nmprep$OCC[nmprep$MDV == 0])
	levels(fitdata$OCCf) <- paste("Occasion",levels(fitdata$OCCf))

	fitdata$GRPf <- factor(fitdata$GRP)
	levels(fitdata$GRPf) <- paste("Group",levels(fitdata$GRPf))

	fitdata$DOSELVLf <- factor(fitdata$DOSELVL)
	levels(fitdata$DOSELVLf) <- paste(c(2.5,2.5,2.5,5,7.5,15,20,25,30,35,50,75),
																		c("mg QD","mg-5mg QD","mg-7.5mg QD",rep("mg QD",9)),sep="")
	fitdata$GENDf <- factor(fitdata$GEND)
	levels(fitdata$GENDf) <- c("F","M")

	fitdata$RACEf <- factor(nmprep$RACE[nmprep$MDV == 0])
	levels(fitdata$RACEf) <- c("Caucasian","Non-Caucasian")

	fitdata$DXCATf <- factor(fitdata$DXCAT)
	levels(fitdata$DXCATf) <- c("CLL","AML","ALL","MM")

	fitdata$CRCL <- nmprep$CRCL2[nmprep$MDV == 0]

#--------------------------------------------------------------------------------------------------
#Diagnostic plots
filename <- "diagnostic_dashboard.png"
png(filename, width = 720, height = 920)

# 4 ggplot2 graphs in a grid layout
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,4)))

  #Plot 1
  max.OBS1 <- max(c(fitdata$PRED, fitdata$DV),na.rm=T)

  plotobj1 <- NULL
  plotobj1 <-  ggplot(fitdata[fitdata$DV>0.001,])
  plotobj1 <- plotobj1 + geom_point(aes(x=PRED, y=DV, colour=GRPf), shape=1)
  plotobj1 <- plotobj1 + geom_abline(aes(x=PRED, y=DV), intercept=0, slope=1, colour="black") #Add line of identity
  plotobj1 <- plotobj1 + geom_smooth(aes(x=PRED, y=DV), method=loess, se=T, colour="red")        #Add loess smoothing line
  plotobj1 <- plotobj1 + scale_x_continuous(name="Population Predicted conc (ug/mL)", limits=c(-0.1,max.OBS1))
  plotobj1 <- plotobj1 + scale_y_continuous(name="Observed conc (ug/mL)", limits=c(-0.1,max.OBS1))
  plotobj1 <- plotobj1 + scale_colour_brewer(name="Dose Level", palette="Set1")
  plotobj1 <- plotobj1 + theme(legend.position="none")
  print(plotobj1, vp=vplayout(1:2,1:2))


  #Plot 2
  max.OBS2 <- max(c(fitdata$IPRED, fitdata$DV),na.rm=T)

  plotobj2 <- NULL
  plotobj2 <-  ggplot(fitdata)
  plotobj2 <- plotobj2 + geom_point(aes(x=IPRED, y=DV, colour=GRPf), shape=1)
  plotobj2 <- plotobj2 + geom_abline(aes(x=IPRED, y=DV), intercept=0, slope=1, colour="black") #Add line of identity
  plotobj2 <- plotobj2 + geom_smooth(aes(x=IPRED, y=DV), method=loess, se=T, colour="red")        #Add loess smoothing line
	plotobj2 <- plotobj2 + geom_abline(aes(x=IPRED, y=DV), intercept=0, slope=2, colour = "darkgreen", linetype = "dashed")
	plotobj2 <- plotobj2 + geom_abline(aes(x=IPRED, y=DV), intercept=0, slope=0.5, colour = "darkgreen", linetype = "dashed")
  plotobj2 <- plotobj2+ scale_x_continuous(name="IPRED conc (ug/mL)", limits=c(-0.1,max.OBS2))
  plotobj2 <- plotobj2+ scale_y_continuous(name="Observed conc (ug/mL)", limits=c(-0.1,max.OBS2))
  plotobj2  <-  plotobj2 + scale_colour_brewer(name="Dose Level", palette="Set1")
  plotobj2  <-  plotobj2 + theme(legend.position="none")
  print(plotobj2, vp=vplayout(1:2,3:4))


  #Plot 3  CWRES vs TAFDE
  max.CWRES <- max(abs(fitdata$CWRES),na.rm=T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1

  plotobj3 <- NULL
  plotobj3 <-  ggplot(fitdata)
  plotobj3 <- plotobj3 + geom_point(aes(x=TAD, y=CWRES, colour=GRPf), shape=1)
  plotobj3 <- plotobj3 + geom_abline(aes(x=TAD, y=CWRES),intercept=0, slope=0, colour="black")  #Add zero line
  plotobj3 <- plotobj3 + geom_smooth(aes(x=TAD, y=CWRES), method=loess, se=T, colour="red")        #Add loess smoothing line
  plotobj3 <- plotobj3+ scale_x_continuous(name="Time after dose (h)")
  plotobj3 <- plotobj3+ scale_y_continuous(name="CWRES", limits=c(-max.CWRES ,max.CWRES))
  plotobj3  <-  plotobj3 + scale_colour_brewer(name="Dose Level", palette="Set1")
  print(plotobj3, vp=vplayout(3,1:4))


  #Plot 4   CWRES vs PRED
  max.CWRES <- max(abs(fitdata$CWRES),na.rm=T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1

  plotobj4 <- NULL
  plotobj4 <-  ggplot(fitdata)
  plotobj4 <- plotobj4 + geom_point(aes(x=PRED, y=CWRES, colour=GRPf), shape=1)
  plotobj4 <- plotobj4 + geom_abline(aes(x=PRED, y=CWRES),intercept=0, slope=0, colour="black")  #Add zero line
  plotobj4 <- plotobj4 + geom_smooth(aes(x=PRED, y=CWRES), method=loess, se=T, colour="red")        #Add loess smoothing line
  plotobj4 <- plotobj4 + scale_x_continuous(name="Population Predicted conc (ug/mL)")
  plotobj4 <- plotobj4 + scale_y_continuous(name="CWRES", limits=c(-max.CWRES ,max.CWRES))
  plotobj4 <-  plotobj4 + scale_colour_brewer(name="Dose Level", palette="Set1")
  print(plotobj4, vp=vplayout(4,1:4))

dev.off()

#Individual diagnostic plots
  to.png.sqr(plotobj1, "DV_PRED")
  to.png.sqr(plotobj2, "DV_IPRED")
  to.png.sqr(plotobj3, "CWRES_TAFDE")
  to.png.sqr(plotobj4, "CWRES_PRED")

#--------------------------------------------------------------------------------------------------
diag.plot<-function(df1,icol,dlow,dup){
filename <- paste("diagnostic_dashboard_",icol,dlow,ifelse(dlow!=dup,dup,""),".png",sep="")
col.interest <- which(names(fitdata)%in%icol)
fitdata <- df1[df1[col.interest]<=dup&df1[col.interest]>=dlow,]
png(filename, width = 720, height = 920)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,4)))

  #Plot 1
	max.OBS1 <- max(c(fitdata$PRED, fitdata$DV),na.rm=T)

  plotobj1 <- NULL
  plotobj1 <-  ggplot(fitdata[fitdata$DV>0.001,])
  plotobj1 <- plotobj1 + geom_point(aes(x=PRED, y=DV), shape=1)
  plotobj1 <- plotobj1 + geom_abline(aes(x=PRED, y=DV), intercept=0, slope=1, colour="black") #Add line of identity
  plotobj1 <- plotobj1 + geom_smooth(aes(x=PRED, y=DV), method=loess, se=T, colour="red")        #Add loess smoothing line
	plotobj1 <- plotobj1 + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
  plotobj1 <- plotobj1+ scale_x_continuous(name="Population Predicted conc (ug/mL)", limits=c(-0.1,max.OBS1))
  plotobj1 <- plotobj1+ scale_y_continuous(name="Observed conc (ug/mL)", limits=c(-0.1,max.OBS1))
  print(plotobj1, vp=vplayout(1:2,1:2))


  #Plot 2
	max.OBS2 <- max(c(fitdata$IPRED, fitdata$DV),na.rm=T)

  plotobj2 <- NULL
  plotobj2 <-  ggplot(fitdata)
  plotobj2 <- plotobj2 + geom_point(aes(x=IPRED, y=DV), shape=1)
  plotobj2 <- plotobj2 + geom_abline(aes(x=IPRED, y=DV), intercept=0, slope=1, colour="black") #Add line of identity
  plotobj2 <- plotobj2 + geom_smooth(aes(x=IPRED, y=DV), method=loess, se=T, colour="red")        #Add loess smoothing line
	plotobj2 <- plotobj2 + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
  plotobj2 <- plotobj2+ scale_x_continuous(name="IPRED conc (ug/mL)", limits=c(-0.1,max.OBS2))
  plotobj2 <- plotobj2+ scale_y_continuous(name="Observed conc (ug/mL)", limits=c(-0.1,max.OBS2))
  print(plotobj2, vp=vplayout(1:2,3:4))


  #Plot 3  CWRES vs TAFDE
  max.CWRES <- max(abs(fitdata$CWRES),na.rm=T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1

  plotobj3 <- NULL
  plotobj3 <-  ggplot(fitdata)
  plotobj3 <- plotobj3 + geom_point(aes(x=TAD, y=CWRES), shape=1)
  plotobj3 <- plotobj3 + geom_abline(aes(x=TAD, y=CWRES),intercept=0, slope=0, colour="black")  #Add zero line
  plotobj3 <- plotobj3 + geom_smooth(aes(x=TAD, y=CWRES), method=loess, se=T, colour="red")        #Add loess smoothing line
  plotobj3 <- plotobj3+ scale_x_continuous(name="Time after dose (h)")
  plotobj3 <- plotobj3+ scale_y_continuous(name="CWRES", limits=c(-max.CWRES ,max.CWRES))
  print(plotobj3, vp=vplayout(3,1:4))


  #Plot 4   CWRES vs PRED
  max.CWRES <- max(abs(fitdata$CWRES),na.rm=T)
  if (max.CWRES < 0.1) max.CWRES <- 0.1

  plotobj4 <- NULL
  plotobj4 <-  ggplot(fitdata)
  plotobj4 <- plotobj4 + geom_point(aes(x=PRED, y=CWRES), shape=1)
  plotobj4 <- plotobj4 + geom_abline(aes(x=PRED, y=CWRES),intercept=0, slope=0, colour="black")  #Add zero line
  plotobj4 <- plotobj4 + geom_smooth(aes(x=PRED, y=CWRES), method=loess, se=T, colour="red")        #Add loess smoothing line
  plotobj4 <- plotobj4 + scale_x_continuous(name="Population Predicted conc (ug/mL)")
  plotobj4 <- plotobj4 + scale_y_continuous(name="CWRES", limits=c(-max.CWRES ,max.CWRES))
  print(plotobj4, vp=vplayout(4,1:4))

dev.off()
	}
	diag.plot(fitdata,"DOSELVL",1,4)
	diag.plot(fitdata,"DOSELVL",5,7)
	diag.plot(fitdata,"DOSELVL",8,9)
	diag.plot(fitdata,"DOSELVL",10,11)
	diag.plot(fitdata,"STUDY",5115,5115)
	diag.plot(fitdata,"STUDY",6003,6003)
	diag.plot(fitdata,"STUDY",8056,8056)
	diag.plot(fitdata,"STUDY",10016,10016)
#--------------------------------------------------------------------------------------------------
#Examine Residual distribution
 #QQ plot
  plotobj <- NULL
	plotobj1 <- NULL
  plotobj <-  ggplot(fitdata)
  plotobj <- plotobj + ggtitle("QQ normal plot of residuals\n")
  plotobj <- plotobj + scale_x_continuous(name="Theoretical", limits=c(-3,3)) + scale_y_continuous(name="Observed", limits=c(-3,3))
  plotobj <- plotobj + geom_abline(aes(sample=CWRES), intercept=0, slope=1, colour="black") #Add line of identity
	plotobj1 <- plotobj + stat_qq(aes(sample=CWRES),geom="point")
  #plotobj

  to.png.sqr(plotobj1,"residual_qq")

 #Condition qq plot on group
  plotobj2 <- NULL
  plotobj2 <- plotobj + stat_qq(aes(sample=CWRES, colour=GRPf),geom="point", alpha=0.5)
  plotobj2  <-  plotobj2 + scale_colour_brewer(name="Group", palette="Set1")
  #plotobj1
  to.png.sqr(plotobj2,"residual_qq_group")

#Condition qq plot on group
  plotobj3 <- NULL
  plotobj3 <- plotobj + stat_qq(aes(sample=CWRES, colour=STUDYf),geom="point", alpha=0.5)
	plotobj3  <-  plotobj3 + scale_colour_brewer(name="Study", palette="Set1")
	#plotobj2
	to.png.sqr(plotobj3,"residual_qq_study")


 #Residual density plot (replaces histogram)
  plotobj <- NULL
  plotobj <-  ggplot(fitdata)
  plotobj <- plotobj + ggtitle("Distribution density of residuals\n")
  plotobj <- plotobj + geom_density(aes(x=CWRES, y=..density..))
  plotobj <- plotobj + scale_x_continuous(name="Residual", limits=c(-4,4))
  plotobj <- plotobj + scale_y_continuous(name="distribution density")
  #plotobj

  to.png.sqr(plotobj,"residual_density")


 #Condition density plot on group
  plotobj <- plotobj + geom_density(aes(x=CWRES, y=..density.., colour=GRPf))
  plotobj  <-  plotobj + scale_colour_brewer(name="Group", palette="Set1")
  #plotobj

  to.png.sqr(plotobj,"residual_density_group")

#--------------------------------------------------------------------------------------------------
#Correlation matrix between eta's
#Subset the eta's from file - make sure they have been tabled in the *.fit file

#ETA correlation--------------------------------------------------------------------------------------------------------------------

#Make a dataframe with only 1 line per subject (avoids duplications of ETA as it is the same for all observations in a subject)
fitdataone <- ddply(fitdata, .(ID), oneperID)
#Can't use this as there are several dose levels per patient

eta.cols <-  grep(glob2rx("ETA*"), names(fitdataone))
etabov.cols <- c("ETA4","ETA5","ETA6","ETA7")  #this is hard-coded for BOV
#eta.cols <- c("ETA4","ETA5")  #hard-coded#
eta.cols <- eta.cols[eta.cols %in% etabov.cols==F]
eta.cols <- names(fitdataone)[eta.cols]
etadata <- subset(fitdataone, select=eta.cols)

if (ncol(etadata)>1)  #more than 1 ETA is scatterplot matrix
{
  plotobj <- NULL
  plotobj <-  ggpairs(etadata, title ="Correlation between random effects")
  #plotobj
  to.png(plotobj,"etascatter")
} else               #only 1 ETA  is density plot
{
  plotobj <- NULL
  plotobj  <- ggplot(data=etadata)
  plotobj <- plotobj + geom_density(aes(x=ETA1, y=..density..), colour="black")
  plotobj <- plotobj +  scale_x_continuous("ETA1")
  plotobj <- plotobj +  scale_y_continuous("distribution density")
  to.png(plotobj,"etascatter")
}

#--------------------------------------------------------------------------------------------------
#Look at ETA versus covariate relationships
#Make a dataframe with only 1 line per subject (avoids duplications of ETA as it is the same for all observations in a subject)
#uses above fitdataone <- ddply(fitdata, .(ID), oneperID)

#eta.cols <-  grep(glob2rx("ETA*"), names(fitdataone))
#etabov.cols <- c("ETA1","ETA2","ETA3")  #this is hard-coded for BOV
#eta.cols <- c("ETA4","ETA5")  #hard-coded#
#eta.cols <- eta.cols[eta.cols %in% etabov.cols==F]

#Get the columns with categorical covariates
covcat.cols <- c("STUDYf","DOSELVLf","GENDf","RACEf","DXCATf","OCCf")

#Get the columns with continuous covariates
covcont.cols <- c("AGE","WT","HT","SECR","CRCL","IBW")

#-------------------------------------------------------------------------------------------------------
#Categorical covariates
covcatdf <- expand.grid(eta.cols,covcat.cols,stringsAsFactors = F)
names(covcatdf) <- c("ETAname","covname")
#covcatdf


#Function to count numbers in a boxplot and return as a coardinate,label pair
boxplot.give.n <- function(x)
{
  return(c(y = median(x), label = length(x)))
}

#Function to plot effect of categorical covariate on ETA
ETACovariatePlotCAT <- function(ETAname,covname)
{

  plotobj <- NULL
  ETAtext <- ETAname  #This could be made an argument to the function
  covtext <- covname  #This could be made an argument to the function
  plotobj <-  ggplot(data=fitdataone)
  plotobj <- plotobj + geom_boxplot(aes_string(x = covname, y = ETAname), position=position_dodge(width=0.9))
  plotobj <- plotobj + stat_summary(aes_string(x = covname, y = ETAname), fun.data = boxplot.give.n, geom = "text", size = 6, colour="red")
  plotobj <- plotobj + scale_x_discrete(paste(covtext))
  plotobj <- plotobj + scale_y_continuous(paste(ETAtext))
  #plotobj  <- plotobj + ggtitle("Final PK model\n")  #legend.position="none",
  #plotobj <- plotobj + ggtitle("Base PK model\n")  #legend.position="none",
  plotobj <- plotobj + ggtitle("ETA PLOT\n")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  #plotobj

  png.file.name <- paste(ETAname,"_vs_",covname,sep="")

  to.png.sqr(plotobj,png.file.name)

}

#Apply the plotting function - stand back and watch the magic!
#This runs the ETACovariatePlotCAT plotting function, taking each row of covcatdf as the input to the function
mdply(covcatdf, ETACovariatePlotCAT)

#-------------------------------------------------------------------------------------------------------
#Continuous covariates
#Make a little dataframe with each row being the arguments for our covariate plotting functions
covcontdf <- expand.grid(eta.cols,covcont.cols,stringsAsFactors = F)
names(covcontdf) <- c("ETAname","covname")
covcontdf

#Function to plot effect of continuous covariate on ETA
ETACovariatePlotCONT <- function(ETAname,covname)
{

  plotobj <- NULL
  ETAtext <- ETAname  #This could be made an argument to the function
  covtext <- covname  #This could be made an argument to the function
  plotobj <-  ggplot(data=fitdataone)
  plotobj <- plotobj + geom_point(aes_string(x = covname, y = ETAname), colour="blue", size=2)
  plotobj <- plotobj + geom_smooth(aes_string(x = covname, y = ETAname), method=loess, se=T, linetype=2)
  plotobj <- plotobj + geom_smooth(aes_string(x = covname, y = ETAname), method=lm, se=F, colour="red")
  plotobj <- plotobj + geom_abline(intercept=0, slope=0, colour="black")
  plotobj <- plotobj + scale_x_continuous(paste(covtext))
  plotobj <- plotobj + scale_y_continuous(paste(ETAtext))
  #plotobj <- plotobj + theme(legend.position="none") + ggtitle("Final PK model\n")
  #plotobj <- plotobj + ggtitle("Base PK model\n")
  plotobj <- plotobj + ggtitle("ETA PLOT\n")
	plotobj <- plotobj + theme(axis.text.x = element_text(angle=45, hjust = 1))
  #plotobj

  png.file.name <- paste(ETAname,"_vs_",covname,sep="")
  to.png.sqr(plotobj,png.file.name)

}

#Apply the plotting function - stand back and watch the magic!
#This runs the ETACovariatePlotCAT plotting function, taking each row of covcontdf as the input to the function
	mdply(covcontdf, ETACovariatePlotCONT)

  param.cols <- c("CL", "V2", "KTR")

  covcatpf <- expand.grid(param.cols,covcat.cols,stringsAsFactors = F)
  names(covcatpf) <- c("ETAname","covname")
	mdply(covcatpf, ETACovariatePlotCAT)

  covcontpf <- expand.grid(param.cols,covcont.cols,stringsAsFactors = F)
  names(covcontpf) <- c("ETAname","covname")
	mdply(covcontpf, ETACovariatePlotCONT)

#Covariate correlation----------------------------------------------------------

	cov.cols <- c(covcat.cols,covcont.cols)
	cov.cols <- c(cov.cols,"CRCL")
	covsubdata <- subset(fitdataone, select=cov.cols)
	plotobj <- NULL
  plotobj <- ggpairs(covsubdata, title ="Correlation between covariates")

	#ggsave("covscatter.png", width=5, height=5, units=c("cm"))

#--------------------------------------------------------------------------------------------------
#Complex ETA grids

	#plotobj1 <- NULL
	#plotobj1 <- ggplot(fitdata)
	#titletext <- expression(atop("ETA1 v's ROUTE",
	#												atop("BY ROUTE",
	#														 "coloured by DOSING")))
	#plotobj1 <- plotobj1 + ggtitle(titletext)
	#plotobj1 <- plotobj1 + geom_point(aes(x=ROUTE, y=ETA1, colour = DOSING), shape=1)
	#plotobj1 <- plotobj1+ scale_x_discrete(name="ROUTE")
	#plotobj1 <- plotobj1+ scale_y_continuous(name="ETA1")
	#plotobj1 <- plotobj1 + facet_wrap(~DOSING)
	#plotobj1

	#png.file.name <- paste("ETA1_grid")
	#to.png.sqr(plotobj1,png.file.name)

	#plotobj2 <- NULL
	#plotobj2 <- ggplot(fitdata)
	#titletext <- expression(atop("ETA2 v's ROUTE",
	#												atop("BY ROUTE",
	#														 "coloured by DOSING")))
	#plotobj2 <- plotobj2 + ggtitle(titletext)
	#plotobj2 <- plotobj2 + geom_point(aes(x=ROUTE, y=ETA2, colour = DOSING), shape=1)
	#plotobj2 <- plotobj2+ scale_x_discrete(name="ROUTE")
	#plotobj2 <- plotobj2+ scale_y_continuous(name="ETA2")
	#plotobj2 <- plotobj2 + facet_wrap(~DOSING)
	#plotobj2

	#png.file.name <- paste("ETA2_grid")
	#to.png.sqr(plotobj2,png.file.name)

#--------------------------------------------------------------------------------------------------
# Plot individual fits and population fits by individual
  plotobj <- NULL
  plotobj <- ggplot(fitdata[fitdata$DV>0.001,])
  titletext <- expression(atop("Individual fits",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobj <- plotobj + ggtitle(titletext)
  plotobj <- plotobj + geom_point(aes(x=TAD, y=DV), shape=1)
  #plotobj <- plotobj + geom_point(aes(x=TIME, y=DV), shape=1)
  plotobj <- plotobj + geom_line(aes(x=TAD, y=PRED), colour = "blue", data=fitdata[fitdata$DV>0.001&fitdata$STUDY==8056,])	#to account for non day 1 dataset
	plotobj <- plotobj + geom_line(aes(x=TAD, y=IPRED), colour = "red", data=fitdata[fitdata$DV>0.001&fitdata$STUDY==8056,])	#to account for non day 1 dataset
	plotobj <- plotobj + geom_line(aes(x=TIME, y=PRED), colour = "blue")		#not TAD due to graphical errors due to two IPREDs
	plotobj <- plotobj + geom_line(aes(x=TIME, y=IPRED), colour = "red")		#not TAD due to graphical errors due to two PREDs
	plotobj <- plotobj + geom_hline(yintercept=0.25, linetype = 2, colour = "darkgreen")
 #plotobj <- plotobj + geom_hline(yintercept=0.3, linetype = 2, colour = "brown")
 #plotobj <- plotobj + geom_hline(yintercept=0.2, linetype = 2, colour = "purple")
	plotobj <- plotobj + facet_wrap(~ ID)
  plotobj <- plotobj+ scale_x_continuous(name="Time after dose (h)", lim = c(0,6))
  plotobj <- plotobj+ scale_y_continuous(name="Lenalidomide (ug/mL)")
  #plotobj

	ggsave("Individual_concs_by_dose_ID_DOSELVL.png", width=90, height=50, units=c("cm"))

	plotobj <- plotobj+ scale_y_log10(name="Lenalidomide Conc (ug/mL)", lim=c(0.01, 5))

	ggsave("Individual_concs_by_dose_ID_DOSELVL_log.png", width=90, height=50, units=c("cm"))

#--------------------------------------------------------------------------------------------------
# Plot 20 best/worst individual fits based on the median prediction errors.
# Not sure I like the use of median PE.  Should really be the sum of the PE's

#First calculate the Prediction Error
  fitdata$PE <- (fitdata$DV-fitdata$IPRED)/fitdata$IPRED
  #head(fitdata)

#Calculate and sort by Median Absolute PE
  #First calculate and order the fitdata by MDAPE
  #Subset for subjects with more than 1 observation to stop graphs crashing
  median.abs <- function(x) median(abs(x), na.rm=T)
  medianPE.df <- summaryBy(PE ~ ID, data=fitdata, FUN=c(median.abs,length))

  medianPE.df.sorted <- orderBy(~ PE.median.abs, data=medianPE.df)
  medianPE.df.sorted <- subset(medianPE.df.sorted, PE.length > 1)

	write.csv(medianPE.df.sorted,file="median.df.sorted.csv", row.names=F)

#Calculate and sort by SSE (sum of PE's squared)
	fitdata$SE <-  fitdata$PE*fitdata$PE
	#fitdata
	sum.abs <- function(x) sum(abs(x), na.rm=T)
  sumSE.df <- summaryBy(SE ~ ID, data=fitdata, FUN=c(sum.abs,length))

  sumSE.df.sorted <- orderBy(~ SE.sum.abs, data=sumSE.df)
  sumSE.df.sorted <- subset(sumSE.df.sorted, SE.length > 1)

	write.csv(sumSE.df.sorted,file="sum.df.sorted.csv", row.names=F)

#Calculate and sort by MSE
	MSE.df.sorted <- sumSE.df.sorted
	MSE.df.sorted <- rename(MSE.df.sorted, c("SE.sum.abs"="PEsum","SE.length"="nData")) #just makes it clearer to me
	MSE.df.sorted$nPar <- 3 # hard coded for the number of fixed effects
	MSE.df.sorted$DegFree <- MSE.df.sorted$nData - MSE.df.sorted$nPar
	MSE.df.sorted$MSE <- MSE.df.sorted$PEsum/MSE.df.sorted$DegFree

  MSE.df.sorted <- orderBy(~ MSE, data=MSE.df.sorted)
  MSE.df.sorted <- subset(MSE.df.sorted, nData > 1)

	write.csv(MSE.df.sorted,file="MSE.df.sorted.csv", row.names=F)

#Plot top 2 fits medianPE---------------------------
  IDtop <- head(subset(medianPE.df.sorted,PE.length>3)$ID, 2)
  fitdatatop <- fitdata[fitdata$ID %in% IDtop,]

	plotobjtop <- NULL
  plotobjtop <- ggplot(fitdatatop)
	titletext <- expression(atop("Best individual fits by medianPE",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobjtop <- plotobjtop + ggtitle(titletext)
	plotobjtop <- plotobjtop + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjtop <- plotobjtop + geom_line(aes(x=TIME, y=PRED), colour = "blue")
	plotobjtop <- plotobjtop + geom_line(aes(x=TIME, y=IPRED), colour = "red")
	plotobjtop <- plotobjtop + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjtop <- plotobjtop + scale_x_continuous(name="Time after dose (h)",lim=c(0,24))
  plotobjtop <- plotobjtop + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjtop <- plotobjtop + facet_grid(~ ID)
	#plotobjtop

  to.png(plotobjtop,"CONC_vs_TAD_by_ID_best_medianPE")

	plotobjtop <- plotobjtop + scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjtop,"LOG_CONC_vs_TAD_by_ID_best_medianPE")

 #Plot bottom 2 fits medianPE------------------------
  IDbottom <- tail(medianPE.df.sorted$ID, 2)
  fitdatabottom <- fitdata[fitdata$ID %in% IDbottom,]

	plotobjbottom <- NULL
  plotobjbottom <- ggplot(fitdatabottom)
	titletext <- expression(atop("Worst individual fits by medianPE",
													atop("Black = observed",
															 "Red  = individual prediction, blue = population prediction")))
	plotobjbottom <- plotobjbottom + ggtitle(titletext)

	plotobjbottom <- plotobjbottom + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjbottom <- plotobjbottom + geom_line(aes(x=TIME, y=PRED), colour = "blue")
	plotobjbottom <- plotobjbottom + geom_line(aes(x=TIME, y=IPRED), colour = "red")
  plotobjbottom <- plotobjbottom + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjbottom <- plotobjbottom + scale_x_continuous(name="Time after dose (h)",lim=c(0,24))
  plotobjbottom <- plotobjbottom + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjbottom <- plotobjbottom + facet_grid( ~ ID)
	#plotobjbottom

  to.png(plotobjbottom,"CONC_vs_TAD_by_ID_worst_medianPE")

	plotobjbottom <- plotobjbottom+ scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjbottom,"LOG_CONC_vs_TAD_by_ID_worst_medianPE")

#Plot top 2 fits MSE---------------------------
  IDtop <- head(subset(MSE.df.sorted,nData>3)$ID, 2)
  fitdatatop <- fitdata[fitdata$ID %in% IDtop,]

	plotobjtop <- NULL
  plotobjtop <- ggplot(fitdatatop)
	titletext <- expression(atop("Best individual fits by MSE",
													atop("Black = observed",
															 "Red = individual prediction, blue = population prediction")))
	plotobjtop <- plotobjtop + ggtitle(titletext)
	plotobjtop <- plotobjtop + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjtop <- plotobjtop + geom_line(aes(x=TAD, y=PRED), colour = "blue")
	plotobjtop <- plotobjtop + geom_line(aes(x=TAD, y=IPRED), colour = "red")
	plotobjtop <- plotobjtop + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
	plotobjtop <- plotobjtop + scale_x_continuous(name="Time after dose (h)")
  plotobjtop <- plotobjtop + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjtop <- plotobjtop + facet_wrap(~ID)
	#plotobjtop

  to.png(plotobjtop,"CONC_vs_TAD_by_ID_best_MSE")

	plotobjtop <- plotobjtop + scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjtop,"LOG_CONC_vs_TAD_by_ID_best_MSE")


#Plot bottom 2 fits MSE------------------------
  IDbottom <- tail(MSE.df.sorted$ID, 2)
  fitdatabottom <- fitdata[fitdata$ID %in% IDbottom,]

	plotobjbottom <- NULL
  plotobjbottom <- ggplot(fitdatabottom)
	titletext <- expression(atop("Worst individual fits by MSE",
													atop("Black = observed",
															 "Red  = individual prediction, blue = population prediction")))
	plotobjbottom <- plotobjbottom + ggtitle(titletext)

	plotobjbottom <- plotobjbottom + geom_point(aes(x=TAD, y=DV), colour = "black")
  plotobjbottom <- plotobjbottom + geom_line(aes(x=TAD, y=PRED), colour = "blue")
	plotobjbottom <- plotobjbottom + geom_line(aes(x=TAD, y=IPRED), colour = "red")
  plotobjbottom <- plotobjbottom + geom_hline(yintercept=0.5, linetype = 2, colour = "darkgreen")
  plotobjbottom <- plotobjbottom + scale_x_continuous(name="Time after dose (h)")
  plotobjbottom <- plotobjbottom + scale_y_continuous(name="Lenalidomide Conc (ug/mL)")
	plotobjbottom <- plotobjbottom + facet_wrap(~ID)
	#plotobjbottom

  to.png(plotobjbottom,"CONC_vs_TAD_by_ID_worst_MSE")

	plotobjbottom <- plotobjbottom+ scale_y_log10(name="Lenalidomide Conc (ug/mL)")

	to.png(plotobjbottom,"LOG_CONC_vs_TAD_by_ID_worst_MSE")
