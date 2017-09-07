# SERVER - Example of using the plot brush with interactive plots
# -----------------------------------------------------------------------------

shinyServer(function(input, output, session) {

# Set up reactive values to keep track of which subset you are up to
  rv <- reactiveValues(
    n = 1,
    flag = Fdata
  )  # reactiveValues
  
# Use this to make a title that informs you of what subset will be used
  output$title <- renderText({
    paste("ID =", sub.id[rv$n])
  })

# Create a subset of the data depending on how many times the user has pressed
# sub.id is located in global.R and is the unique values for the column you
# plan to subset your data by
  Rdata <- reactive({
    fitdata[fitdata$ID == sub.id[rv$n], ]
  })  # Rdata()

# Plot the subsetted data
# The plot is just a regular plot unless made interactive in ui.R
  output$plot <- renderPlot({
    data <- Rdata()[-which(Rdata()$MDV == 1), ]
    plotobj <- ggplot(data, aes(x = TAD, y = DV)) +
      geom_point() +
      theme_bw() +
      scale_x_continuous("Time after dose (hrs)", lim = c(0,24)) +
      scale_y_log10("DV (ug/mL)")
    plotobj
  })  # plot
  
# Create reactive object containing the brushed data
  Bdata <- reactive({
    brushedPoints(Rdata(), input$plot_brush)
  })

# Output the info gained from using brush on plot points
  output$info <- renderPrint({
    Bdata()[sub.verb]
  })  # output$info

# Setup the previous next output to cycle through plots
  output$nextPrev <- renderUI({
    if (rv$n == length(sub.id)) {
      prevFun(1)
    }  # if on final tab
    else if (rv$n == 1) {
      nextFun(10)  # 10 so that it is created in the same spot
    }  # if on first tab
    else {
      div(prevFun(1), nextFun(8))
    }  # otherwise
  })  # renderUI

# Use the left and right buttons to cycle through the plot subsets
# Do this by using obsereEvents that change what part of the subset is next
  observeEvent(input$prevTab, {
    rv$n <- rv$n - 1
  })  # observeEvent prevTab

  observeEvent(input$nextTab, {
    rv$n <- rv$n + 1
  })  # observeEvent nextTab
  
# Observe for when user flags the rows and then save this data
  observeEvent(input$flag, {
    rv$flag[Bdata()$X] <- 1
  })
  
# Save copy of nmprep with flagged rows
  Sdata <- reactive({
    data <- nmprep
    data$FLAG <- rv$flag
    return(data)
  })
  
# Output the flagged rows to keep a track of them
  output$flagged <- renderPrint({
    nmprep[which(rv$flag == 1), c("#ID", "TIME", "TAD", "STUDY")]
  })
  
# Observe for when user wants to save the work they have done
  observeEvent(input$save, {
    write.csv(Sdata(), "E:/Hughes/Git/LEN_PK/select_app/Data/nmprep_flagged.csv",
      quote = F, row.names = F)
  })

# When window containing the app is closed, stop the app
  session$onSessionEnded(function() {
    stopApp()
  })  # onSessionEnded
  
# Debug Console
  observe(label = "console", {
    if(input$console != 0) {
      options(browserNLdisabled = TRUE)
      saved_console <- ".RDuetConsole"
      if (file.exists(saved_console)) load(saved_console)
      isolate(browser())
      save(file = saved_console, list = ls(environment()))
    }
  })  # observe
})
