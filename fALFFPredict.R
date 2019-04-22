
# fALFFPredict.R

library(caret)

source('~/R/Stroke Aphasia/Treatment Response Prediction/import.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/impute.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/predict.R')

falff.predict <- function(site, num.imputations){
  
  # Import data 
  data <- import()
  data <- data[-c(40,47),] # Remove dropped subjects
  data <- data[data$Site == site,]
  data <- data[,c(70:90)]
  
  # Remove subjects with missing scans
  subjects.missing <- which(is.na(data$FALFF1Pre))
  data <- data[-subjects.missing,]
  
  # Split data into predictors, target
  falff <- data[,2:21]
  deficit.post <- data[,1]
  
  # Elastic net regression
  r.val <- score.predict.imp(falff, deficit.post, num.imputations, site)

  return(r.val)  
  
}