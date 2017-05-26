# server.R script for PMB_lenPopPK_app
# ------------------------------------------------------------------------------
shinyServer(function(input, output, session) {

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
			CRCL = as.numeric(input$crcl),
			FFM = as.numeric(input$ffm)
  	)

  	oral.dose.times <- 0
  	oral.dose.data <- input.conc.data[input.conc.data$time %in% oral.dose.times, ]
  	oral.dose.data$amt <- input$dose
  	oral.dose.data$evid <- 1
  	oral.dose.data$rate <- 0
  	oral.dose.data$cmt <- 1

  	input.conc.data <- rbind(input.conc.data, oral.dose.data)
  	input.conc.data <- input.conc.data[with(input.conc.data, order(input.conc.data$ID, input.conc.data$time)), ]

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
		if (input$ci90) {
	    plotobj1 <- plotobj1 + stat_summary(aes(x = time, y = IPRE), data = Rconc(),
				geom = "ribbon", fun.ymin = CI90lo, fun.ymax = CI90hi,
				fill = cPalette[rv$n + 1], alpha = 0.2
			)
		}

		if (length(rv$Sconc) != 0) {
			plotobj1 <- plotobj1 + stat_summary(
				aes(x = TIME, y = IPRE, colour = palette), data = rv$Sconc,
				geom = "line", fun.y = median, size = 1)
      if (input$ci90) {
				plotobj1 <- plotobj1 + stat_summary(
					aes(x = TIME, y = IPRE, fill = palette), data = rv$Sconc,
					geom = "ribbon", fun.ymin = CI90lo, fun.ymax = CI90hi, alpha = 0.2
				)
			}
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

	output$clValue <- renderValueBox({
		valueBox(
			paste(signif(median(Rconc()$CL), 3), "L/h"),
			tags$p("Clearance", style = "font-size: 20px"),
			icon = icon("log-out", lib = "glyphicon"),
			color = "purple"
		)  # valueBox.cl
	})

	output$vdValue <- renderValueBox({
		valueBox(
			paste(signif(median(Rconc()$V1), 3), "L"),
			tags$p("Volume of Distribution", style = "font-size: 20px"),
			icon = icon("refresh"),
			color = "green"
		)  # valueBox.vd
	})

	# Close the R session when browser closes
	session$onSessionEnded(function(){
	 stopApp()
	})  #endsession
})
