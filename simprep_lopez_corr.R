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

  sex_prob <- 10/15  # for males

  dose_cat <- c(25, 15, 10, 5)
  dose_prob <- c(6/15, 5/15, 2/15, 2/15)

# -----------------------------------------------------------------------------
# Create resampling for normal distribution if outside the desired range
  trunc_rnorm <- function(n, mean, sd, range, log = F) {
    lower <- min(range)
    upper <- max(range)
    # Convert mean and sd if necessary
    if (log) {  # transform mu and sigma to work in log domain
      mu <- log(mean^2/sqrt((sd^2 + mean^2)))
      sigma <- sqrt(log(sd^2/mean^2 + 1))
      # https://au.mathworks.com/help/stats/lognstat.html
    } else {
      mu <- mean
      sigma <- sd
    }
    # Sample from normal distribution
    out <- rnorm(n, mu, sigma)
    if (log) { out <- exp(out) }
    # Determine values to be resampled
    repeat {
      resample <- out < lower | out > upper
      if (!any(resample)) break  # if none to resample exit loop
      replace <- which(resample)
      new <- rnorm(length(replace), mu, sigma)
      if (log) { new <- exp(new) }
      out[replace] <- new
    }
    out
  }

  trunc_mvrnorm <- function(n, mean, sd, corr_mat, lower, upper, log = F) {
    require(MASS)
    require(MBESS)
    mat_dim <- length(mean)  # frequently used term
    # Check for errors
    if (length(sd) != mat_dim | length(lower) != mat_dim | length(upper) != mat_dim) {
      stop("mean, sd, lower and upper all must have the same length")
    } else if (!is.matrix(corr_mat)) {
      stop("correlation matrix must be of class matrix")
    } else if (any(!dim(corr_mat) %in% mat_dim)) {
      stop("correlation matrix dimensions must match length of mean")
    }
    # Convert mean and sd if necessary
    if (log) {  # transform mu and sigma to work in log domain
      mu <- log(mean^2/sqrt((sd^2 + mean^2)))
      sigma <- sqrt(log(sd^2/mean^2 + 1))
      # https://au.mathworks.com/help/stats/lognstat.html
    } else {
      mu <- mean
      sigma <- sd
    }
    # Sample from multivariate normal distribution
    cov_mat <- cor2cov(corr_mat, sigma)
    if (log) {  # take normally distributed error terms for transformation
      mvr_mat <- mvrnorm(n, double(mat_dim), cov_mat)
      out <- apply(mvr_mat, 2, exp) %*% diag(mean)
    } else {  # take the normal distribution
      out <- mvrnorm(n, mean, cov_mat)
    }
    # Determine values to be resampled
    repeat {
      # mapply - check if values are outside of their respective range
      # matrix - to coerce from vector to matrix
      # apply - to see which rows have any values that need resampling
      resample <- apply(MARGIN = 1, FUN = any,
        matrix(ncol = mat_dim,
          mapply(function(val, low, upp) {
            val < low | val > upp
          }, out, rep(lower, each = n), rep(upper, each = n))
        )  # matrix
      )  # apply
      if (!any(resample)) break  # if none to resample exit loop
      if (log) {
        mvr_mat <- matrix(ncol = mat_dim,
          mvrnorm(length(which(resample)), double(mat_dim), cov_mat)
        )  # matrix used as length(which(resample)) == 1 returns a vector
        out[resample, ] <- apply(mvr_mat, 2, exp) %*% diag(mean)
      } else {
        out[resample, ] <- mvrnorm(length(which(resample)), mean, cov_mat)
      }
    }
    out
  }

  # Testing trunc_mvrnorm
  # corr <- matrix(c(
  #   1.0, 0.7, 0.5,
  #   0.7, 1.0, 0.0,
  #   0.5, 0.0, 1.0), 3, 3)
  # test <- trunc_mvrnorm(100, c(74.6, 165, 1.8), c(18.2, 8, 0.2), corr,
  # c(51, 154, 1.5), c(109, 182, 2.3), F)
  # plot(test[,2] ~ test[,3])

# -----------------------------------------------------------------------------
# Simulate patient demographics
# Set number of simulated patients and simulation times
  nstudy <- 1000
  nid_perstudy <- 15
  nid <- nstudy*nid_perstudy
  id_seq <- 1:nid
  times <- unique(c(seq(0, 4, 0.25), seq(4, 8, 0.5), seq(8, 24, 1)))

# Sample doses and sex given demographic proportions
  dose <- sample(dose_cat, nid, replace = T, prob = dose_prob)
  sex <- rbinom(nid, 1, sex_prob)

# Set up correlation matrix for height and weight
# Age does not require correlation as Cockroft-Gault handles this
  r11 <- 1  # correlation between weight and weight
  r12 <- 0.7  # correlation between weight and height
  r22 <- 1  # correlation between height and height
  corr <- matrix(c(
    r11, r12,
    r12, r22
  ), 2, 2)  # symmetrical correlation matrix

# Simulate covariate values
  corrsim <- trunc_mvrnorm(nid, c(wt_mean, ht_mean), c(wt_sd, ht_sd), corr,
    c(wt_range[1], ht_range[1]), c(wt_range[2], ht_range[2]), log = T)
  wt <- corrsim[, 1]
  ht <- corrsim[, 2]
  age <- trunc_rnorm(nid, age_mean, age_sd, age_range, log = F)
  secr <- trunc_rnorm(nid, secr_mean, secr_sd, secr_range, log = T)
  crcl_sex <- sex*1.23
  crcl_sex[crcl_sex == 0] <- 1.04
  crcl <- (140 - age)*wt*crcl_sex/secr

# Resample until crcl is within the range outlined in the paper
  repeat {
    resample <- crcl < min(crcl_range) | crcl > max(crcl_range)
    if (!any(resample)) break  # if none to resample exit loop
    replace <- which(resample)
    new_nid <- length(replace)
    new_corrsim <- trunc_mvrnorm(new_nid, c(wt_mean, ht_mean), c(wt_sd, ht_sd),
      corr, c(wt_range[1], ht_range[1]), c(wt_range[2], ht_range[2]), log = T)
    wt[replace] <- new_corrsim[, 1]
    ht[replace] <- new_corrsim[, 2]
    age[replace] <- trunc_rnorm(new_nid, age_mean, age_sd, age_range, log = F)
    secr[replace] <- trunc_rnorm(new_nid, secr_mean, secr_sd, secr_range, log = T)
    crcl[replace] <- (140 - age[replace])*wt[replace]*crcl_sex[replace]/secr[replace]
  }

# Perform transformations on covariate values to attain BSA and FFM
  bsa <- 0.007184*wt^0.425*ht^0.725
  ffm_c1 <- sex*6.68e3
  ffm_c1[ffm_c1 == 0] <- 8.78e3
  ffm_c2 <- sex*216
  ffm_c2[ffm_c2 == 0] <- 244
  bmi <- wt/(ht/100)^2
  ffm <- 9.27e3 * wt / (ffm_c1 + ffm_c2 * bmi)

# -----------------------------------------------------------------------------
# Prepare simulation dataset
  simprep_amt <- data.frame(
    ID = id_seq,
    STUDY = ceiling(id_seq/15),
    AMT = dose,
    SEX = sex,
    WT = round(wt, 2),
    HT = round(ht, 2),
    CRCL2 = round(crcl, 2),
    BSA = round(bsa, 2),
    FFM = round(ffm, 2)
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
  filename_out <- paste(output_dir,"simdata_lopez_mvr15000.csv",sep="/")
  write.csv(simprep, file = filename_out, quote = F, row.names = F)
