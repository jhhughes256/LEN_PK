### Trial using Box-Cox and Manly transforms
# -----------------------------------------------------------------------------
# Set up environment
  rm(list = ls(all = T))

# Load libraries
  library(ggplot2)
  library(plyr)

# Create log-normal distribution
  var <- 0.36
  eta <- rnorm(1000, 0, sqrt(var))
  
# Apply Box-Cox transform
  par_bx <- 0.1
  phi <- exp(eta)
  eta_bx <- (phi**par_bx-1)/par_bx
  
# Apply Manly transform
  par_man <- 1.5
  eta_man <- (phi*par_man-1)/par_man
  
# Compare the distributions
  eta_df <- data.frame(
    type = rep(c("norm", "boxcox", "manly"), each = 1000),
    eta = c(phi, exp(eta_bx), exp(eta_man))
  )
  
  p <- NULL
  p <- ggplot(data = eta_df)
  p <- p + geom_density(aes(x = eta))
  p <- p + scale_x_continuous(lim = c(0, 5))
  p <- p + facet_wrap(~type, ncol = 3)
  p

  