library(shiny)
library(ggplot2)

fitdata <- read.csv("E:/Hughes/Data/PK/REDO/RUN063_DES_1C8TAwAP_PPV_CORCLVKA_FFM.nm7/run063_des_1c8tawap_ppv_corclvka_ffm.fit.csv")

sub.id <- sort(unique(fitdata$X.ID))

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
