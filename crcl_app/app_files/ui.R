# Description
# ------------------------------------------------------------------------------
  fixedPage(
    h3("Lenalidomide Population Pharmacokinetic Model"),
    selectInput("mod",
      "Choose Model",
      choices = list(
        "Structural Model" = 1,
        "Covariate Model" = 2
      ),
      selected = 1
    ),  # selectInput
    hr(),
    plotOutput("concPlot",width = 600),  # Concentration-time profile output
    fixedRow(
      column(5,
        selectInput("crcl",
          "Creatinine Clearance (ml/min)",
          choices = list(10, 30, 60, 90),
          selected = 90
        )  # sliderInput
      ),  # column
      column(5,
        numericInput("dose",
          "Dose (mg):",
          value = 25,
          step = 5
        )  # numericInput 
      ),  # column
      column(2,
        actionButton("save", "Save"),  # actionButton
        actionButton("clear", "Clear")  # actionButton
      ) # column
    ), # fixedRow
    fixedRow(
      div("Press the save button to keep the current curve. Then proceed to
      change the covariates to compare the saved curve with the reactive curve."),
      strong("Saving more than 8 plots is not recommended as I don't know
      more than 8 colours."),
      hr(),
      checkboxInput("log", "Plot concentrations on a log scale")
    ),  # fixedRow
    align = "center"
  )  # fixedPage
