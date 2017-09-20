# Simulate correlated weights and heights
# ------------------------------------------------------------------------------
# Load required packages
  library(MASS)	# mvrnorm function
  library(MBESS)	# cor2cov function

# ------------------------------------------------------------------------------
# Number of individuals to be simulated
  nid <- 1000
  ID.seq <- 1:nid	# Create a sequence of ID numbers

# ------------------------------------------------------------------------------
# Assign mean and standard deviation values to weight and height
# Weight
  mean.WT <- 70	# kg
  sd.WT <- 0.2
# Height
  mean.HT <- 170	# cm
  sd.HT <- 0.1

# Set up correlation matrix for weight and height
  r11 <- 1	# Correlation between weight and weight
  r12 <- 0.7	# Correlation between weight and height
  r22 <- 1	# Correlation between height and height
  corr <- matrix(c(
    r11,r12,
    r12,r22
  ),2,2)	# Correlation matrix needs to be symmetrical

# Convert to correlation matrix to a variance-covariance matrix
  cov.mat <- cor2cov(cor.mat = corr,sd = c(sd.WT,sd.HT))

# ------------------------------------------------------------------------------
# Simulate values for weight and height using the variance-covariance matrix
# mvrnorm is a multivariate normal distribution random number generator
# i.e., it is not for a log normal distribution
# The easiest way would be to generate normally distributed "error" terms
# And then transform them like we would log-normally distributed ETAs
  sim.values <- mvrnorm(n = nid,mu = c(0,0),Sigma = cov.mat)
  df <- data.frame(ID = ID.seq,
    WT = mean.WT*exp(sim.values[,1]),	# WT "error" terms were the first column of sim.values matrix hence [,1]
    HT = mean.HT*exp(sim.values[,2]))

# Plot distribution of weight
  plot(hist(df$WT))

# Plot distribution of height
  plot(hist(df$HT))

# Plot weight values versus height values
  plot(df$WT ~ df$HT)

# Perform linear regression to check the correlation
  lm.result <- lm(df$WT ~ df$HT)
  summary(lm.result)
  # Adjusted R-squared = 0.4826
  # r = sqrt(0.4826) = 0.7 = r12 from above
