# Description
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(dplyr)	# New plyr - required for mrgsolve
	library(mrgsolve)	# Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# Set colour palette
	cPalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
		"#FFFF33", "#A65628", "#F781BF")

# ------------------------------------------------------------------------------
# Set number of individuals that make up the 95% prediction intervals
	n <- 1000
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x, probs = 0.025)
	CI95hi <- function(x) quantile(x, probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x, probs = 0.05)
	CI90hi <- function(x) quantile(x, probs = 0.95)
# Set seed for reproducible numbers
	set.seed(123456)

	TIME <- seq(from = 0,to = 24,by = 0.25)

# Source the models
	source("model_base.R")
	source("model_cov.R")

# ------------------------------------------------------------------------------
  shinyServer(function(input, output) {

		Rconc <- reactive({
			ID <- 1:n
			ID2 <- sort(c(rep(ID, times = length(TIME))))
			time <- rep(TIME, times = length(ID))
    	input.conc.data <- data.frame(
    		ID = ID2,
    		time,
    		amt = 0,
    		evid = 0,
    		rate = 0,
    		cmt = 1,
				CRCL = as.numeric(input$crcl)
    	)

    	oral.dose.times <- 0
    	oral.dose.data <- input.conc.data[input.conc.data$time %in% oral.dose.times, ]
    	oral.dose.data$amt <- input$dose
    	oral.dose.data$evid <- 1
    	oral.dose.data$rate <- 0
    	oral.dose.data$cmt <- 1

    	input.conc.data <- rbind(input.conc.data, oral.dose.data)
    	input.conc.data <- input.conc.data[with(input.conc.data, order(input.conc.data$ID, input.conc.data$time)), ]

			if (input$mod == 1) {
				mod <- mod.base
			} else if (input$mod == 2) {
				mod <- mod.cov
			}

    	conc.data <- mod %>%
        data_set(input.conc.data) %>%
        mrgsim()
    	conc.data <- as.data.frame(conc.data)
    })  # Rconc

		rv <- reactiveValues(
      n = 0,  # additional plots
      Sconc = data.frame(Cp = NULL, palette = NULL)  # saved concentrations
    )  # rv

		observeEvent(input$save, {
      rv$n <- rv$n + 1
      rv$Sconc <- rbind(
        rv$Sconc,
        data.frame(
          IPRE = Rconc()$IPRE,
					TIME = Rconc()$time,
          palette = factor(rv$n)
        )
      )
    })  #observeEvent

		observeEvent(input$clear, {
      rv$n <- 0
      rv$Sconc <- data.frame(Cp = NULL, palette = NULL)
    })  #observeEvent

    output$concPlot <- renderPlot({
      plotobj1 <- NULL
      plotobj1 <- ggplot()
      plotobj1 <- plotobj1 + stat_summary(aes(x = time, y = IPRE), data = Rconc(),
				geom = "line", fun.y = median, colour = cPalette[rv$n + 1], size = 1)
      plotobj1 <- plotobj1 + stat_summary(aes(x = time, y = IPRE), data = Rconc(),
				geom = "ribbon", fun.ymin = CI90lo, fun.ymax = CI90hi,
				fill = cPalette[rv$n + 1], alpha = 0.2)

			if (length(rv$Sconc) != 0) {
				plotobj1 <- plotobj1 + stat_summary(
					aes(x = TIME, y = IPRE, colour = palette), data = rv$Sconc,
					geom = "line", fun.y = median, size = 1)
	      plotobj1 <- plotobj1 + stat_summary(
					aes(x = TIME, y = IPRE, fill = palette), data = rv$Sconc,
					geom = "ribbon", fun.ymin = CI90lo, fun.ymax = CI90hi, alpha = 0.2)
			}

			plotobj1 <- plotobj1 + scale_colour_manual(values = cPalette)
			plotobj1 <- plotobj1 + scale_fill_manual(values = cPalette)
			plotobj1 <- plotobj1 + theme(legend.position="none")

      plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)", lim = c(0,24))
			if (input$log) {
				plotobj1 <- plotobj1 + scale_y_log10("Concentration (mg/L)\n")
			} else {
				plotobj1 <- plotobj1 + scale_y_continuous("Concentration (mg/L)\n")
			}
      print(plotobj1)
    })  # renderPlot
  })
