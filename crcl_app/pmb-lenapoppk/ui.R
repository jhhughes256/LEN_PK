# ui.R script for PMB_lenPopPK_app
# ------------------------------------------------------------------------------
  header <- dashboardHeader(
    titleWidth = 400,
    title = "Lenalidomide Model Demonstration"
  )  # dashboardHeader

  sidebar <- dashboardSidebar(
    disable = TRUE
  )  # dashboardSidebar

  body <- dashboardBody(
    # Concentration-time plot
    fluidRow(
      box(width = 12, align = "center",
        h2("Mean Concentration Time Profile from 500 Simulations"),
        plotOutput("concPlot"),
        fluidRow(style = "font-size: 18px",
          checkboxInput("ci90",
            "Plot upper and lower 90% confidence intervals"
          ),
          checkboxInput("log",
            "Plot concentrations on a log scale"
          )
        )  # fluidRow.plotOptions
      )  # box.plot
    ),  # fluidRow.plot
    fluidRow(
      # Inputs
      box(width = 8, style = "padding-left: 20px",
        # Column 1
        column(width = 10, style = "font-size: 18px",
          fluidRow(
            selectInput("ffm",
              "Fat Free Mass (kg):",
              choices = list(45, 55, 65),
              selected = 55
            )  # numericInput.ffm
          ),  # fluidRow.ffm
          fluidRow(
            selectInput("crcl",
              "Creatinine Clearance (ml/min):",
              choices = list(30, 60, 90, 120),
              selected = 90
            )  # sliderInput.crcl
          ),  # fluidRow.crcl
          fluidRow(
            numericInput("dose",
              "Dose (mg):",
              value = 25,
              step = 5
            )  # numericInput.dose
          )  # fluidRow.dose
        ),  # column.input
        column(width = 2, align = "center",
          style = "padding: 40px",
          fluidRow(style = "padding: 10px",
            actionButton(
              "save", "", icon("floppy-o"),
              style = "
                color: #fff;
                background-color: #337ab7;
                border-color: #2e6da4;
                font-size: 25px;
                padding: 5px 15px 8px
              "  # style.saveButton
            )  # actionButton.save
          ),  # fluidRow.save
          fluidRow(style = "padding: 10px",
            actionButton("clear", "", icon("trash-o"),
              style = "
                color: #fff;
                background-color: #b83333;
                border-color: #a32e2e;
                font-size: 25px;
                padding: 5px 15px 8px
              "  # style.clearButton
            )  # actionButton.clear
          )  # fluidRow.clear
        )  # column.saveclear
      ),  # box.options
      # Text Outputs
      column(width = 4, style = "padding-top: 10px; font-size: 18px",
        fluidRow(
          valueBoxOutput("clValue", width = 12)
        ),
        fluidRow(
          valueBoxOutput("vdValue", width = 12)
        )
      )
    ),  # fluidRow.options
    tags$head(
      tags$style(
        HTML('
          .skin-blue .main-header .logo {
            background-color: #3c8dbc;
          }
          .skin-blue .main-header .logo:hover {
            background-color: #3c8dbc;
          }
        ')
      )
    )
  )  # dashboardBody

  dashboardPage(
    header, sidebar, body
  )
