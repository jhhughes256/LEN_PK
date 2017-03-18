# SERVER - Example of using the plot brush with interactive plots
# -----------------------------------------------------------------------------

shinyServer(function(input, output, session) {

# Set up reactive values to keep track of which subset you are up to
  rv <- reactiveValues(n = 1)

# Create a subset of the data depending on how many times the user has pressed
# sub.cyl is located in global.R and is the unique values for the column you
# plan to subset your data by
  Rdata <- reactive(
    mtcars[mtcars$cyl == sub.cyl[rv$n], ]
  )  # Rdata()

# Plot the subsetted data
# The plot is just a regular plot unless made interactive in ui.R
  output$plot <- renderPlot({
    plotobj <- ggplot(Rdata(), aes(x = wt, y = mpg)) +
      geom_point() +
      theme_bw()
    plotobj
  })  # plot

# Output the info gained from using brush on plot points
  output$info <- renderPrint({
    brushedPoints(mtcars, input$plot_brush)
  })  # output$info

# Setup the previous next output to cycle through plots
  output$nextPrev <- renderUI({
    if (rv$n == length(sub.cyl)) {
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

# When window containing the app is closed, stop the app
  session$onSessionEnded(function() {
    stopApp()
  })  # onSessionEnded
})
