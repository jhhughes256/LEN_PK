# Prepare dataset for simulation
# -----------------------------------------------------------------------------
# Dataset will then be simulated using my model and the Guglieri-Lopez model
# -----------------------------------------------------------------------------
# Prepare work environment
# Remove any previous objects in the workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set the working directory
  master_dir <- "E:/Hughes/Data"
  scriptname <- "sim_prep"
  setwd(master_dir)

# Load required packages
  library(MASS)  # mvrnorm
  library(MBESS)  # cor2cov

# Organise working and output directories
  working_dir <- paste(master_dir,"RAW_Clinical",sep="/")
  workspacefilename <- paste(getwd(),"/",scriptname,".RData", sep="")

  output_dir <- paste(working_dir,"/",scriptname,"_Output",sep="")
  if(!file.exists(output_dir)){
	  dir.create(output_dir)
  }

# -----------------------------------------------------------------------------
# Patient demographics for Guglieri-Lopez dataset
  age_mean <- 63.8
  age_sd <- 8.5
  age_range <- c(53, 84)

  wt_mean <- 74.6
  wt_sd <- 18.2
  wt_range <- c(51, 109)

  ht_mean <- 165
  ht_sd <- 8
  ht_range <- c(154, 182)

  secr_mean <- 0.91*88.42
  secr_sd <- 0.27*88.42
  secr_range <- c(0.7*88.42, 1.84*88.42)

  bsa_mean <- 1.8
  bsa_sd <- 0.2
  bsa_range <- c(1.5, 2.3)

  crcl_mean <- 89.4
  crcl_sd <- 38.1
  crcl_range <- c(31.8, 175.2)

# -----------------------------------------------------------------------------
# Create resampling for normal distribution if outside the desired range
  trunc_rnorm <- function(n, mean, sd, range, log = F) {
    lower <- min(range)
    upper <- max(range)
    if (log) {
      mu <- log(mean^2/sqrt((sd^2 + mean^2)))
      sigma <- sqrt(log(sd^2/mean^2 + 1))
    } else {
      mu <- mean
      sigma <- sd
    }
    out <- rnorm(n, mu, sigma)
    if (log) { out <- exp(out) }
    repeat {
      resample <- out < lower | out > upper
      if (!any(resample)) break
      replace <- which(resample)
      new <- rnorm(length(replace), mu, sigma)
      if (log) { new <- exp(new) }
      out[replace] <- new
    }
    out
  }

# -----------------------------------------------------------------------------
# Simulate patient demographics
# Set number of simulated patients and simulation times
  nid <- 1000
  id_seq <- 1:nid
  times <- unique(c(seq(0, 4, 0.25), seq(4, 8, 0.5), seq(8, 24, 1)))

# Set doses and proportion of all doses
  dose <- c(
    rep(25, round(6/15*nid)),
    rep(15, round(5/15*nid)),
    rep(10, round(2/15*nid)),
    rep(5, round(2/15*nid))
  )
  if (length(dose) > nid) {
    excess <- 1:(length(dose) - nid) + nid
    dose <- dose[-excess]
  } else if (length(dose) < nid) {
    missing <- (nid + 1) - 1:(nid - length(dose))
    dose[missing] <- unique(dose)[1:length(missing)]
  }

# Simulate covariate values
  crcl <- trunc_rnorm(nid, crcl_mean, crcl_sd, crcl_range, log = T)
  crcl <- round(crcl, 2)

  bsa <- trunc_rnorm(nid, bsa_mean, bsa_sd, bsa_range, log = T)
  bsa <- round(bsa, 2)

# -----------------------------------------------------------------------------
# Prepare simulation dataset
  simprep_amt <- data.frame(
    ID = id_seq,
    AMT = dose,
    BSA = bsa,
    CRCL = crcl
  )
  simprep_dv <- data.frame(
    ID = rep(id_seq, each = length(times)),
    TIME = times,
    EVID = 0,
    DV = 1,
    CMT = 2,
    MDV = 0
  )
  simprep <- merge(simprep_dv, simprep_amt)

# Fix up columns
  simprep$EVID[simprep$TIME == 0] <- 1
  simprep$DV[simprep$TIME == 0] <- "."
  simprep$CMT[simprep$TIME == 0] <- 1
  simprep$MDV[simprep$TIME == 0] <- 1
  simprep$AMT[simprep$TIME != 0] <- "."

# Save as .csv to directory
  names(simprep)[1] <- "#ID"
  filename_out <- paste(output_dir,"simdata_lopez.csv",sep="/")
  write.csv(simprep, file = filename_out, quote = F, row.names = F)
