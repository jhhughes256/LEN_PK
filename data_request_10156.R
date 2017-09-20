# Request for creatinine clearance data
# Requires patient numbers and date of corresponding doses
# -----------------------------------------------------------------------------
# Set the working directory
  master.dir <- "E:/Hughes/Data/RAW_Clinical"
  scriptname <- "data_request_10156"

# Load libraries
  library(plyr)
  library(readxl)

# Source utility functions file
  source(paste0("E:/Hughes/functions_utility.r"))

# -----------------------------------------------------------------------------
# Load data
  file.name.in1 <- paste0(master.dir,"/rawdata-lena_10156/pk data with new patient data.xlsx")
  datanew.pkd <- as.data.frame(read_excel(file.name.in1, sheet=1))  #pk data

  file.name.in2 <- paste0(master.dir,"/rawdata-lena_10156/Data_PK analysis.xlsx")
  datapk.demog <- as.data.frame(read_excel(file.name.in2, sheet=1))  #demographic data
  datapk.dose <- as.data.frame(read_excel(file.name.in2, sheet=6))  #actual dose times

# Adjust column names for usability in R
  names(datanew.pkd) <- str_replace_all(names(datanew.pkd),"[ ()#]",".")
  head(datanew.pkd)

# -----------------------------------------------------------------------------
# Check for errors in patient numbers
  with(datanew.pkd, table(Patient.Number, useNA = "always"))
  # There appears to be one patient number entered as NA. There is also a
  # 8834-23 and 8834-23b, these are not the same patient.
  # 8834-23b is thought to be 8834-29

  na.row <- with(datanew.pkd, which(is.na(Patient.Number)))
  datanew.pkd[(na.row-2):(na.row+2),]
  # The NA row is the final row can be removed

  datanew.pkd <- datanew.pkd[-dim(datanew.pkd)[1], ]

# Extract patient numbers and sample date
  dataone <- ddply(datanew.pkd, .(Patient.Number), function(x) { x[1, ] })
  dataout <- dataone[c(2, 3, 4, 5)]

# Save file
  filename.out <- paste(master.dir, "datacheck_clin_10156_Output", "missing_crcl_10156.csv", sep = "/")
  write.csv(dataout, file=filename.out, row.names=FALSE)
