
library(missForest)

# Imputes missing data using missForest package 
impute.rf <- function(data){
  
  # Convert input data (likely data frame) into a matrix
  data <- data.matrix(data)
  
  # Impute the data with missForest package (random forests approach) (non-verbose)
  invisible(capture.output(data.imputed <- missForest(data, variablewise = TRUE, mtry = 67)))

  # Return imputed data
  return(data.imputed$ximp)
  
}