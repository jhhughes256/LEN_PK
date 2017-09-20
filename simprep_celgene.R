# Prepare dataset for simulation
# -----------------------------------------------------------------------------
# Dataset will then be simulated using my model and the Celgene model
# May be difficult to simulate the dataset used in this model...
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

# Create function for determining mean from median and range
  median2mean <- function(n, median, range) {
    lower <- min(range)
    upper <- max(range)
    if (nid <= 25) {
      mean <- (lower + 2*median + upper)/4
    } else {
      mean <- median
    }
    if (nid <= 15) {
      sd <- sqrt(((lower - 2*median + upper)/4 + (upper - lower)^2)/12)
    } else if (nid > 70) {
      sd <- (upper - lower)/4
    } else {
      sd <- (upper - lower)/6
    }
    c(mean = mean, sd = sd)
  }

# -----------------------------------------------------------------------------
# Patient demographics for Celgene datasets
# All median + range combos converted as per
# Hozo SP et. al 2005
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-5-13
# Richardson PG et. al 2002 (MM-001 & MM-002)  n = 34
# Missing 14 patients
  rich_nid <- 34

  rich_age_median <- 59
  rich_age_range <- c(40, 69)
  rich_age_meansd <- median2mean(rich_nid, rich_age_median, rich_age_range)

  rich_wt_median <- 82
  rich_wt_range <- c(50, 118)
  rich_wt_meansd <- median2mean(rich_nid, rich_wt_median, rich_wt_range)

  rich_crcl_median <- 101
  rich_crcl_range <- c(65, 155)
  rich_crcl_meansd <- median2mean(rich_nid, rich_crcl_median, rich_crcl_range)

# Iida S et. al 2010 (MM-017)  n = 12
# Missing 2 patients
  iida_nid <- 12

  iida_age_median <- 63
  iida_age_range <- c(43, 66)
  iida_age_meansd <- median2mean(iida_nid, iida_age_median, iida_age_range)

  iida_wt_median <- 59
  iida_wt_range <- c(48, 75)
  iida_wt_meansd <- median2mean(iida_nid, iida_wt_median, iida_wt_range)

  iida_crcl_median <- 91
  iida_crcl_range <- c(63, 135)
  iida_crcl_meansd <- median2mean(iida_nid, iida_crcl_median, iida_crcl_range)

# Hou J et. al 2012 (MM-021)  n = 9
# Missing one patient that didn't have crcl > 60ml/min
  hou_nid <- 9

  hou_age_median <- 55
  hou_age_range <- c(44, 68)
  hou_age_meansd <- median2mean(hou_nid, hou_age_median, hou_age_range)

  hou_wt_median <- 65
  hou_wt_range <- c(54, 84)
  hou_wt_meansd <- median2mean(hou_nid, hou_wt_median, hou_wt_range)

  hou_crcl_median <- 95
  hou_crcl_range <- c(63, 154)
  hou_crcl_meansd <- median2mean(hou_nid, hou_crcl_median, hou_crcl_range)

#


# -----------------------------------------------------------------------------
# Simulate patient demographics
# Set number of simulated patients and simulation times
  nid <- 1000
  id_seq <- 1:nid
  times <- unique(c(seq(0, 4, 0.25), seq(4, 8, 0.5), seq(8, 24, 1)))

# Set doses and proportion of all doses
  dose <- c(
    rep(25, round(6/83*nid)),
    rep(15, round(5/83*nid)),
    rep(10, round(2/83*nid)),
    rep(5, round(2/83*nid))
  )
  if (length(dose) > nid) {
    excess <- 1:(length(dose) - nid) + 1000
    dose <- dose[-excess]
  } else if (length(dose) < nid) {
    missing <- 1001 - 1:(nid - length(dose))
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
  filename_out <- paste(output_dir,"simdata_celgene.csv",sep="/")
  write.csv(simprep, file = filename_out, quote = F, row.names = F)
