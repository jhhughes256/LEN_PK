# R script for simulating a population from a generic 2-compartment model
# ------------------------------------------------------------------------------
# Load package libraries
	library(ggplot2)	# Plotting
	library(grid)	# Plotting
	library(dplyr)	# New plyr - required for mrgsolve
	library(mrgsolve)	# Metrum differential equation solver for pharmacometrics
# Define a custom ggplot2 theme
	theme_bw2 <- theme_set(theme_bw(base_size = 16))

# ------------------------------------------------------------------------------
# Set number of individuals that make up the 95% prediction intervals
	n <- 1000
# 95% prediction interval functions - calculate the 2.5th and 97.5th percentiles
	CI95lo <- function(x) quantile(x,probs = 0.025)
	CI95hi <- function(x) quantile(x,probs = 0.975)
# 90% prediction interval functions - calculate the 5th and 95th percentiles
	CI90lo <- function(x) quantile(x,probs = 0.05)
	CI90hi <- function(x) quantile(x,probs = 0.95)
# Set seed for reproducible numbers
	set.seed(123456)

	TIME <- seq(from = 0,to = 24,by = 0.25)

# ------------------------------------------------------------------------------
# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# This compiled model is used for simulating n individuals and their concentration-time profiles
	code <- '
  	$INIT // Initial conditions for compartments
  		DEPOT = 0,  // Depot - dose enters the system here
			TRAN1 = 0,  // Transit compartment 1
			TRAN2 = 0,  // Transit compartment 2
			TRAN3 = 0,  // Transit compartment 3
			TRAN4 = 0,  // Transit compartment 4
			TRAN5 = 0,  // Transit compartment 5
			TRAN6 = 0,  // Transit compartment 6
			TRAN7 = 0,  // Transit compartment 7
  		CENT = 0,  // Central
  		AUC = 0  // Area under the curve compartment

  	$PARAM  // Population parameters
  		POPCL = 10,  // Clearance, L/h
  		POPV1 = 72.1,  // Volume of central compartment, L
  		POPKTR = 13.3,  // Absorption rate constant, h^-1

		  // Default covariate values for simulation
			STUDY = 0,  // Patient study
			WT = 70,  // Total body weight (kg)
			HT = 170,  // Patient height (cm)
			SEX = 1  // Patient sex

  	$OMEGA  // Omega covariance block
      block = TRUE
  		labels = s(BSV_CL,BSV_V1)
  		0.4651  // BSV for CL
  		0.2480 0.2970  // BSV for V1

		$OMEGA  // Omega variance
			labels = s(BSV_KTR)
  		0.3684  // BSV for KTR

  	$SIGMA  // Sigma
      block = FALSE
  		labels = s(ERR_ME,ERR_LO,ERR_HI)
			0.2510  // Proportional error combined
			0.1163  // Proportional error 5115 6003
  		0.3831  // Proportional error 8156 10016

  	$MAIN  // Determine covariate values
			double BMI = WT/pow(HT*0.01,2);
		  double FFM = 9270*WT/(6680+216*BMI);
			if(SEX == 0) FFM = 9270*WT/(8780+244*BMI);
			// Individual parameter values
  		double CL = POPCL*pow(FFM/55,0.75)*exp(BSV_CL);
  		double V1 = POPV1*(FFM/55)*exp(BSV_V1);
  		double KTR = POPKTR*exp(BSV_KTR);
			// Proportional residual unexplained variability
			double ERR_PRO = ERR_ME;
			if(STUDY >= 5115) ERR_PRO = ERR_LO;
			if(STUDY >= 8156) ERR_PRO = ERR_HI;

  	$ODE  // Differential equations
  		dxdt_DEPOT = -KTR*DEPOT;
			dxdt_TRAN1 = KTR*DEPOT -KTR*TRAN1;
			dxdt_TRAN2 = KTR*TRAN1 -KTR*TRAN2;
			dxdt_TRAN3 = KTR*TRAN2 -KTR*TRAN3;
			dxdt_TRAN4 = KTR*TRAN3 -KTR*TRAN4;
			dxdt_TRAN5 = KTR*TRAN4 -KTR*TRAN5;
			dxdt_TRAN6 = KTR*TRAN5 -KTR*TRAN6;
			dxdt_TRAN7 = KTR*TRAN6 -KTR*TRAN7;
  		dxdt_CENT = KTR*TRAN7 -CL/V1*CENT;
  		dxdt_AUC = CENT/V1;

  	$TABLE  // Determine individual predictions
      double IPRE = CENT/V1;
  		double DV = IPRE*(1+ERR_PRO);

  	$CAPTURE  // Capture output
      IPRE DV CL V1 KTR
	'
	# Compile the model code
	mod <- mcode("popDEMO",code)
	# There is opportunity to simply update model parameters after the model code has been compiled

# ------------------------------------------------------------------------------
# Simulate concentration-time profiles for the population
	ID <- 1:n
	ID2 <- sort(c(rep(ID,times = length(TIME))))
	time <- rep(TIME,times = length(ID))
	input.conc.data <- data.frame(
		ID = ID2,
		time,
		amt = 0,
		evid = 0,
		rate = 0,
		cmt = 1
	)

	oral.dose.times <- 0
	oral.dose.data <- input.conc.data[input.conc.data$time %in% oral.dose.times,]
	oral.dose.data$amt <- 25
	oral.dose.data$evid <- 1
	oral.dose.data$rate <- 0
	oral.dose.data$cmt <- 1

	input.conc.data <- rbind(input.conc.data,oral.dose.data)
	input.conc.data <- input.conc.data[with(input.conc.data, order(input.conc.data$ID,input.conc.data$time)),]

	conc.data <- mod %>% data_set(input.conc.data) %>% mrgsim()
	conc.data <- as.data.frame(conc.data)	# Convert to a data frame so that it is more useful for me!

# ------------------------------------------------------------------------------
# Plot results
	plotobj1 <- NULL
	plotobj1 <- ggplot(conc.data)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "line",fun.y = median,colour = "red",size = 1)
	plotobj1 <- plotobj1 + stat_summary(aes(x = time,y = IPRE),geom = "ribbon",fun.ymin = "CI90lo",fun.ymax = "CI90hi",fill = "red",alpha = 0.3)
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)",lim = c(0,24))
	# plotobj1 <- plotobj1 + scale_y_log10("Concentration (mg/L)\n")
	plotobj1 <- plotobj1 + scale_y_continuous("Concentration (mg/L)\n")
	print(plotobj1)
