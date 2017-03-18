# UI - Example of using the plot brush with interactive plots
# -----------------------------------------------------------------------------

fluidPage(
  h2("Plot Brush Example"),
  plotOutput("plot",
    brush = "plot_brush",
    height = 250, width = 600
  ),  # plotOutput
  fluidRow(
    uiOutput("nextPrev")
  ),  # fluidRow
  hr(),
  fluidRow(
    verbatimTextOutput("info")
  ),  # fluidRow
  align = "center"
)
