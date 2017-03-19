# UI - Example of using the plot brush with interactive plots
# -----------------------------------------------------------------------------

fluidPage(
  textOutput("title", h2),
  actionButton("console","Debug Console"),
  plotOutput("plot",
    brush = brushOpts("plot_brush", resetOnNew = T),
    height = 250, width = 600
  ),  # plotOutput
  fluidRow(
    uiOutput("nextPrev")
  ),  # fluidRow
  hr(),
  fluidRow(
    verbatimTextOutput("info")
  ),  # fluidRow
  fluidRow(
    actionButton("flag",
      "Flag Row(s)"
    )  # actionButton
  ),  # fluidRow
  hr(),
  fluidRow(
    verbatimTextOutput("flagged")
  ),
  fluidRow(
    actionButton("save",
      "Save"
    )  # actionButton
  ),
  align = "center"
)
