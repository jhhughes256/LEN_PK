# Define the model parameters and equations
	# Using mrgsolve - analytical solutions
	# This compiled model is used for simulating n individuals and their concentration-time profiles
	# Cannot have whitespace before $ blocks
  # Simulation uses RUN063_DES_1C8TAwAP_PPV_CORCLVKA_FFM
# ------------------------------------------------------------------------------

code <- '
$INIT // Initial conditions for compartments
  DEPOT = 0,  // Depot - dose enters the system here
  TRAN1 = 0,  // Transit compartment 1
  TRAN2 = 0,  // Transit compartment 2
  TRAN3 = 0,  // Transit compartment 3
  CENT = 0,  // Central
  AUC = 0  // Area under the curve compartment

$PARAM  // Population parameters
  POPCL = 14,  // Clearance, L/h
  POPV1 = 57,  // Volume of central compartment, L
  POPKTR = 6.55,  // Absorption rate constant, h^-1

  // Covariate effects
  COV1 = 0.0065,	// Effect of creatinine clearance on clearance
  COV2 = 1.37,  // Effect of body surface area on volume

  // Default covariate values for simulation
  STUDY = 0,  // Patient study
	BSA = 1.73,  // Fat free mass (kg)
  CRCL = 90  // Creatinine clearance (umol/L)

$OMEGA  // Omega variance
  labels = s(BSV_CL,BSV_V1,BSV_KTR)
  0.0961  // BSV for CL
  0.04  // BSV for V1
  0.3721  // BSV for KTR

$SIGMA  // Sigma
  block = FALSE
  labels = s(ERR_POI)
  0.0144  // Poisson error combined

$MAIN  // Individual parameter values
  double CL = POPCL*(1+COV1*(CRCL-90))*exp(BSV_CL);
  double V1 = POPV1*(1+COV2*(BSA-1.73))*exp(BSV_V1);
  double KTR = POPKTR*exp(BSV_KTR);

$ODE  // Differential equations
  dxdt_DEPOT = -KTR*DEPOT;
  dxdt_TRAN1 = KTR*DEPOT -KTR*TRAN1;
  dxdt_TRAN2 = KTR*TRAN1 -KTR*TRAN2;
  dxdt_TRAN3 = KTR*TRAN2 -KTR*TRAN3;
  dxdt_CENT = KTR*TRAN3 -CL/V1*CENT;
  dxdt_AUC = CENT/V1;

$TABLE  // Determine individual predictions
  double IPRE = CENT/V1;
  double DV = IPRE+pow(IPRE,0.5)*ERR_POI;

$CAPTURE  // Capture output
  IPRE DV CL V1 KTR
'

# Compile the model code
mod.lopez <- mcode("popLOPEZ",code)
