
library(caret)
library(dplyr)
library(ggplot2)

source('~/R/Stroke Aphasia/Treatment Response Prediction/impute.R')

# Elastic net linear regression: output = f(input)
score.predict <- function(input, output){
  
  # Enable leave-one-out cross-validation, bounds on outcome prediction
  control <- trainControl(method="LOOCV", predictionBounds = c(TRUE, TRUE))
  
  # Set range of values for hyperparameters
  fraction.vals <- ((0:100)/100) ** 2  # Ratio of L1 to L2 regularization
  lambda.vals <- 2^(16:-16) # Total regularization penalty
  
  # Make a grid of hyperparameter combinations
  hp.grid <- expand.grid(fraction = fraction.vals, lambda = lambda.vals)
  
  # Convert inputs to data frame
  input <- as.data.frame(input)
  
  # Train an elastic net regression model
  model <- train(input, output, method = "enet", 
                 trControl = control, tuneGrid = hp.grid)
  
  # Pick hyperparameters to optimize correlation, otherwise RMSE
  model <- choose.params(model)
  
  return(model)
  
}

# Update optimal hyperparameters based on correlation
choose.params <- function(model){
  
  # For each row in results
  results <- model$results
  best.cor <- 0
  
  for(idx in 1:nrow(results)){
    
    # Get the hyperparameter combination
    fraction <- model$results$fraction[idx]
    lambda <- model$results$lambda[idx]
    
    # Get the predictions from that combination
    predictions <- model$pred[(model$pred$lambda == lambda & model$pred$fraction == fraction),1:2]
    
    # Calculate the correlation between predicted/actual values
    correlation <- cor(predictions)[2]
    
    # If it is larger than the current best correlation
    if(correlation > best.cor){
      
      # Update the model parameter tuning
      model$bestTune$fraction <- fraction
      model$bestTune$lambda <- lambda
      best.cor <- correlation
      
    }
  }
  
  # Return the model with updated hyperparameters
  return(model)
    
}

# Bootstrap a two-sided confidence interval on model correlation
bootstrap.ci <- function(model, percentile, iterations){
  
  # Get optimal hyperparameters from model
  fraction <- model$bestTune$fraction
  lambda <- model$bestTune$lambda
  
  # Get the model results using the optimal hyperparameters
  results <- model$pred[(model$pred$lambda == lambda & model$pred$fraction == fraction),1:2]
  
  # Bootstrap the predictions
  correlations <- integer(iterations)
  for(i in 1:iterations){
    
    # Resample the prediction results
    results.resampled <- sample_n(results, size = nrow(results), replace = TRUE)
    
    # Save the resampled correlation
    correlations[i] <- cor(results.resampled)[2]
  }
  
  # Sort the bootstrapped correlations
  correlations <- sort(correlations)

  # Calculate median and confidence bounds
  lb.index <- round(iterations * (1 - percentile) / 2)
  ub.index <- iterations - lb.index
  
  lb <- correlations[lb.index]
  ub <- correlations[ub.index]
  med <- median(correlations)
  
  # Return the median and confidence bounds
  return(c(lb, med, ub))
  
}

# Create a performance plot comparing actual vs. predicted values
plot.predict <- function(actual, predicted, title, r.val, mad.val){
  
  # Convert R and N values to character strings
  r.label = paste('R =', as.character(round(r.val, digits = 3)))
  mad.label = paste('MAD =', as.character(round(mad.val, digits = 2)))
  
  # Calculate display limits
  xmin <- min(min(predicted,actual))
  xmax <- max(max(predicted,actual))
  xlabel <- xmin + 0.2 * (xmax - xmin)
  ylabel <- xmax - 0.1 * (xmax - xmin)
  
  # Plot
  p <- ggplot(mapping = aes(x = predicted, y = actual)) +
    theme_classic(base_size = 40) +
    geom_point(size = 10) +
    geom_abline(intercept = 0, slope = 1, size = 3, linetype = "dashed") +
    geom_smooth(size = 5, color = "black", se = FALSE, linetype = "solid") +
    ggtitle(title) +
    annotate("text", x = xlabel, y = xmax, label = r.label, size = 13) +
    #annotate("text", x = xlabel, y = ylabel, label = mad.label, size = 13) +
    labs(x = "Predicted Score", y = "Actual Score") + 
    xlim(xmin,xmax) + ylim(xmin,xmax) + coord_fixed()
  
  print(p)
  
}

# Assess convergence of correlation across imputations
var.converge <- function(var.vals, title, ylabel){
  
  # Set up analysis
  num.imputations <- length(var.vals)
  running.med <- integer(num.imputations)
  
  # Calculate a running median
  for(i in 1:num.imputations){
    running.med[i] <- median(var.vals[1:i])
  }
  
  # Plot
  p <- ggplot(mapping = aes(x = 1:num.imputations, y = running.med)) +
    theme_classic(base_size = 40) +
    geom_line(size = 3) +
    labs(x = "Number of Imputations", y = ylabel) + 
    ggtitle(title)

  print(p)
  
}

# Assess model performance with missing data
score.predict.imp <- function(input, output, num.imp, site){
  
  # Store imputation results
  output.imp <- array(dim = c(num.imp, length(output)))
  
  # Store r interval values
  r.ubs <- integer(num.imp)
  r.lbs <- integer(num.imp)
  r.meds <- integer(num.imp)
  mads <- integer(num.imp)
  
  # Store prediction results
  pred.imp <- array(dim = c(num.imp, length(output)))
  
  # For each imputation
  for(imp in 1:num.imp){
    
    # Impute the missing data 
    data <- data.frame(input, output)
    data <- impute.rf(data)
    
    # Store imputed values
    output.imp[imp,] <- data[,ncol(data)]
    
    # Elastic net regression
    model <- score.predict(data[,1:(ncol(data) - 1)], data[,ncol(data)])
    
    # Constrain predictions to output range
    model$pred$pred[model$pred$pred > 1] = 1
    model$pred$pred[model$pred$pred < 0] = 0
    
    # Bootstrap a confidence interval for correlation
    r.vals <- bootstrap.ci(model, 0.95, 1000)
    r.lbs[imp] <- r.vals[1]
    r.meds[imp] <- r.vals[2]
    r.ubs[imp] <- r.vals[3]

    # Get optimal hyperparameters from model
    fraction <- model$bestTune$fraction
    lambda <- model$bestTune$lambda
    
    # Get the model results using the optimal hyperparameters
    results <- model$pred[(model$pred$lambda == lambda & model$pred$fraction == fraction),1:2]
    mads[imp] <- median(abs(results$pred - results$obs))
    pred.imp[imp,] <- results$pred
    
  }
  
  # Summarize performance
  predicted <- apply(pred.imp, 2, median)
  actual <- apply(output.imp, 2, median)
  r.val <- median(r.meds)
  mad.val <- median(mads)

  # Get the title
  if(site == "BU") {title = "Anomia"}
  if(site == "NU") {title = "Agrammatism"}
  if(site == "JHU") {title = "Dysgraphia"}
  
  # Plot correlation convergence
  var.converge(r.meds, title, "Median Correlation")
  
  # Plot model performance
  plot.predict(actual, predicted, title, r.val, mad.val)
  
  # Print summary statistics
  cat("\n", "--- Prediction Summary Statistics ---", "\n")
  cat("Correlation (Lower Bound): ", as.character(round(median(r.lbs), digits = 3)), "\n")
  cat("Correlation (Median): ", as.character(round(median(r.meds), digits = 3)), "\n")
  cat("Correlation Upper Bound: ", as.character(round(median(r.ubs), digits = 3)), "\n")
  cat("Median Absolute Deviation: ", as.character(round(median(mads), digits = 3)), "\n")
  
  return(r.val)

}
