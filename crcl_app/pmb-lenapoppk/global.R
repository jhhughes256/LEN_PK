# global.R script for PMB_lenPopPK_app
# ------------------------------------------------------------------------------
# Load package libraries
	library(shinydashboard)
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(dplyr)	# New plyr - required for mrgsolve
	library(mrgsolve)	# Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 20))

# Set colour palette
	cPalette <- c("#3C8DBC", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00",
		"#FFFF33", "#A65628", "#F781BF")

# ------------------------------------------------------------------------------
# Set number of individuals that make up the 95% prediction intervals
	n <- 500
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x, probs = 0.025)
	CI95hi <- function(x) quantile(x, probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x, probs = 0.05)
	CI90hi <- function(x) quantile(x, probs = 0.95)
# Set seed for reproducible numbers
	set.seed(123456)

	TIME <- seq(from = 0,to = 24,by = 0.1)

# Source the models
	source("model_cov.R")
