
# VariableImportance.R
# Backwards variable selection

library(caret)

source('~/R/Stroke Aphasia/Treatment Response Prediction/import.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/impute.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/predict.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/(fALFF Predict)/fALFFPredict.R')
source('~/R/Stroke Aphasia/Treatment Response Prediction/(Behavioral Predict)/BehaviorPredict.R')

behavioral.vars <- function(site, num.imp){
  
  # Import data 
  data <- import()
  data <- data[-c(40,47),] # Remove dropped subjects
  data <- data[data$Site == site,c(9:13,15:33,35:37,42,70)]
  
  # Partition data into predictors, deficit
  baseline.tests <- data[,1:28]
  deficit.post <- data[,29]
  num.vars <- ncol(baseline.tests)
  
  # Estimate variable importance by backselection
  is.included <- !integer(num.vars)
  variable.order <- integer(num.vars)
  r.meds <- integer(num.vars)
  
  # For each iteration of backwards selection 
  for(iter in num.vars:2){
    
    start.time <- Sys.time()
    cat("Variables Remaining:", iter, "\n")
    
    # Store the median model correlation
    var.r.meds <- integer(num.vars)
    
    # For each remaining variable
    for(var in which(is.included)){
      
      # Cannot remove pre-Tx deficit measure
      if(var == ncol(baseline.tests)){
        var.r.meds[var] <- -1
        next
      }
      
      # Remove it from the model input
      vars.idx <- is.included
      vars.idx[var] <- FALSE
      vars.idx <- which(vars.idx)
      
      # Store correlations across inputations
      var.cors <- integer(num.imp)
      
      # For each imputation
      for(imp in 1:num.imp){
        
        # Impute the missing data 
        data <- data.frame(baseline.tests, deficit.post)
        data <- impute.rf(data)
        
        # If only two variables remain use linear regression
        if(length(which(is.included)) == 2){
          
          control <- trainControl(method="LOOCV", predictionBounds = c(TRUE, TRUE))
          model <- train(as.data.frame(data[,vars.idx]), 
                         data[,ncol(data)], method = "lm", trControl = control)
          var.cors[imp] <- cor(model$pred$pred, model$pred$obs)
          next
          
        }
        
        # Elastic net regression
        model <- score.predict(data[,vars.idx], data[,ncol(data)])

        # Constrain predictions to output range
        model$pred$pred[model$pred$pred > 1] = 1
        model$pred$pred[model$pred$pred < 0] = 0
        
        # Get optimal hyperparameters from model
        fraction <- model$bestTune$fraction
        lambda <- model$bestTune$lambda
        
        # Get the model results using the optimal hyperparameters
        results <- model$pred[(model$pred$lambda == lambda & model$pred$fraction == fraction),1:2]
        pred <- results$pred
        obs <- results$obs
        
        # Calculate and save correlation
        var.cors[imp] <- cor(pred, obs)
        
      }
      
      # Save the median correlation
      var.r.meds[var] <- median(var.cors)
      
    }
    
    # Remove the variable which maximizes the model correlation
    removed.idx <- which.max(var.r.meds)
    is.included[removed.idx] <- FALSE
    variable.order[iter] <- removed.idx
    r.meds[iter] <- var.r.meds[removed.idx]
    
    print(Sys.time() - start.time)
  }
  
  # Calculate the correlation with all variables
  r.vars.all <- behavior.predict(site, num.imp)
  
  # Fill in the last variable
  variable.order[1] <- setdiff(1:num.vars, variable.order)
  
  # Calculate importance
  r2.meds <- r.meds ** 2
  var.imps <- integer(num.vars)
  var.imps[1:(num.vars-1)] <- r2.meds[2:num.vars]
  var.imps[num.vars] <- r.vars.all ** 2
  var.imps <- var.imps - r2.meds
  
  # Reorder and convert to percentages
  var.imps[variable.order] <- var.imps
  var.imps <- var.imps * 100
  var.names <- colnames(baseline.tests)
  
  # Print output
  cat("Variable Importance Results:", site, "\n")
  print(data.frame(var.names, var.imps))
  
  return(var.imps)
  
}

falff.vars <- function(site, num.imp){
  
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
  num.vars <- ncol(falff)

  # Estimate variable importance by backselection
  is.included <- !integer(num.vars)
  variable.order <- integer(num.vars)
  r.meds <- integer(num.vars)
  
  # For each iteration of backwards selection 
  for(iter in num.vars:2){
    
    start.time <- Sys.time()
    cat("Variables Remaining:", iter, "\n")

    # Store the median model correlation
    var.r.meds <- integer(num.vars)

    # For each remaining variable
    for(var in which(is.included)){
      
      # Remove it from the model input
      vars.idx <- is.included
      vars.idx[var] <- FALSE
      vars.idx <- which(vars.idx)

      # Store correlations across inputations
      var.cors <- integer(num.imp)
      
      # For each imputation
      for(imp in 1:num.imp){
        
        # Impute the missing data 
        data <- data.frame(falff, deficit.post)
        data <- impute.rf(data)
        
        # If only two variables remain use linear regression
        if(length(which(is.included)) == 2){
          
          control <- trainControl(method="LOOCV", predictionBounds = c(TRUE, TRUE))
          model <- train(as.data.frame(data[,vars.idx]), 
                         data[,ncol(data)], method = "lm", trControl = control)
          var.cors[imp] <- cor(model$pred$pred, model$pred$obs)
          next
          
        }
        
        # Elastic net regression
        model <- score.predict(data[,vars.idx], data[,ncol(data)])
        
        # Constrain predictions to output range
        model$pred$pred[model$pred$pred > 1] = 1
        model$pred$pred[model$pred$pred < 0] = 0
        
        # Get optimal hyperparameters from model
        fraction <- model$bestTune$fraction
        lambda <- model$bestTune$lambda
        
        # Get the model results using the optimal hyperparameters
        results <- model$pred[(model$pred$lambda == lambda & model$pred$fraction == fraction),1:2]
        pred <- results$pred
        obs <- results$obs
        
        # Calculate and save correlation
        var.cors[imp] <- cor(pred, obs)
        
      }
      
      # Save the median correlation
      var.r.meds[var] <- median(var.cors)
      
    }
    
    # Remove the variable which maximizes the model correlation
    removed.idx <- which.max(var.r.meds)
    is.included[removed.idx] <- FALSE
    variable.order[iter] <- removed.idx
    r.meds[iter] <- var.r.meds[removed.idx]
    
    print(Sys.time() - start.time)
  }
  
  # Calculate the correlation with all variables
  r.vars.all <- falff.predict(site, num.imp)
  
  # Fill in the last variable
  variable.order[1] <- setdiff(1:num.vars, variable.order)
  
  # Calculate importance
  r2.meds <- r.meds ** 2
  var.imps <- integer(num.vars)
  var.imps[1:(num.vars-1)] <- r2.meds[2:num.vars]
  var.imps[num.vars] <- r.vars.all ** 2
  var.imps <- var.imps - r2.meds
  
  # Reorder and convert to percentages
  var.imps[variable.order] <- var.imps
  var.imps <- var.imps * 100
  var.names <- colnames(falff)
  
  # Print output
  cat("Variable Importance Results:", site, "\n")
  print(data.frame(var.names, var.imps))
  
  return(var.imps)

}