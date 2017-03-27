# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# This compiled model is used for simulating n individuals and their concentration-time profiles
	# Cannot have whitespace before $ blocks
# ------------------------------------------------------------------------------

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
	POPCL = 12  // Clearance, L/h
	POPV1 = 68.8,  // Volume of central compartment, L
	POPKTR = 13.5,  // Absorption rate constant, h^-1

	// Covariate effects
	COV1 = 0.224,	// Effect of creatine clearance on clearance

  // Default covariate values for simulation
	STUDY = 0,  // Patient study
	FFM = 55,  // Fat free mass (kg)
  CRCL = 90  // Creatinine clearance (umol/L)

$OMEGA  // Omega covariance block
  block = TRUE
	labels = s(BSV_CL,BSV_V1)
	0.2959  // BSV for CL
	0.2469 0.2841  // BSV for V1

$OMEGA  // Omega variance
	labels = s(BSV_KTR)
	0.3672  // BSV for KTR

$SIGMA  // Sigma
  block = FALSE
	labels = s(ERR_PRO)
	0.1849  // Proportional error combined

$MAIN  // Determine covariate values
	// Individual parameter values
	double CL = POPCL*pow(FFM/55,0.75)*pow(CRCL/90,COV1)*exp(BSV_CL);
	double V1 = POPV1*(FFM/55)*exp(BSV_V1);
	double KTR = POPKTR*exp(BSV_KTR);

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
	mod.cov <- mcode("popCOV",code)
