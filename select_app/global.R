library(shiny)
library(ggplot2)

fitdata <- read.csv("E:/Hughes/Data/PK/REDO/RUN063_DES_1C8TAwAP_PPV_CORCLVKA_FFM.nm7/run063_des_1c8tawap_ppv_corclvka_ffm.fit.csv")

nmprep <- read.csv("E:/Hughes/Data/PK/REDO/COV24/nmprep_allstudies.csv")
names(nmprep)[1] <- "#ID"

Fdata <- double(length(fitdata$ID))
#Fdata[which(fitdata$TIME > 48 & fitdata$STUDY == 6003 & fitdata$MDV == 0)] <- 1

sub.id <- sort(unique(fitdata$ID))

sub.verb <- c("ID", "TIME", "TAD", "DV", "CWRES", "WRES", "CL", "V2", "KTR")

prevFun <- function(x) {
  column(1,
    actionButton("prevTab",
      icon("angle-double-left",
        class = "fa-2x"
      ),
      width = 70
    ),  # actionButton
    offset = x
  )  # column
}  # prevFun

nextFun <- function(x) {
  column(1,
    actionButton("nextTab",
      icon("angle-double-right",
        class = "fa-2x"
      ),
      width = 70
    ),  # actionButton
    offset = x
  )  # column
}  # nextFun
