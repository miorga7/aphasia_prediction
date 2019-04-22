
# BehaviorPredict.R

source('~/R/Stroke Aphasia/Treatment Response Prediction/import.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/impute.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/predict.R')

behavior.predict <- function(site, num.imputations){
  
  # Import data 
  data <- import()
  data <- data[-c(40,47),] # Remove dropped subjects
  data <- data[data$Site == site,c(9:13,15:33,35:37,42,70)]
  
  # Partition data into predictors, deficit
  baseline.tests <- data[,1:28]
  posttx.deficit <- data[,29]
  
  # Elastic net regression with missing data imputation
  r.val <- score.predict.imp(baseline.tests, posttx.deficit, num.imputations, site)

  return(r.val)
  
}
