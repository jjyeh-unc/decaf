apply_decaf = function(data, classifier){
  
  # "data": a dataframe with unique official gene symbols as rownames
  # "classifier": the DeCAF classifier containing essential objects
  
  ## Extract Gene Games
  genes = rownames(data)
  
  ## Extract Classifier 
  fit = classifier$fit
  if(is.null(fit$beta)) "Classifier Does Not Have Coefficients Assigned to beta"
  
  ## Keep only gene info for genes in classifier 
  data1 = data[genes %in% classifier$TSPs,]
  rnames = rownames(data1)
  data1 = matrix(as.numeric(as.matrix(data1)), ncol = ncol(data1))
  rownames(data1) = rnames
  if(nrow(data1) != length(unique(classifier$TSPs))){
    print(classifier$TSPs[!classifier$TSPs %in% genes])
    stop("genes missing")
  }
  
  ## See which of gene pair has higher expression
  indmat = matrix(-1, ncol(data1), nrow(classifier$TSPs))
  for(i in 1:nrow(classifier$TSPs)){
    p1 = which(rownames(data1) == classifier$TSPs[i,1])
    p2 = which(rownames(data1) == classifier$TSPs[i,2])
    indmat[,i] = (data1[p1,] > data1[p2,])^2
  }
  
  ## Calculate probability of permCAF
  X=cbind(rep(1, nrow(indmat)), indmat)
  trainingPrediction = exp(X%*%c(fit$beta))/(1+exp(X%*%c(fit$beta)))
  
  ## Obtain DeCAF subtype
  classification = c("restCAF","permCAF")[(trainingPrediction >= 0.5)^2 + 1]
  
  ## Obtain Grade of DeCAF Subtype
  guess = rep(1, length(trainingPrediction))
  guess[trainingPrediction < .1] = "Strong restCAF"
  guess[trainingPrediction >= .1 & trainingPrediction < .4] = "Likely restCAF"
  guess[trainingPrediction >= .4 & trainingPrediction < .5] = "Lean restCAF"
  guess[trainingPrediction >= .5 & trainingPrediction < .6] = "Lean permCAF"
  guess[trainingPrediction >= .6 & trainingPrediction < .9] = "Likely permCAF"
  guess[trainingPrediction >= .9 ] = "Strong permCAF"
  
  ## Put results together into dataframe
  final = data.frame(DeCAF_prob = trainingPrediction, DeCAF = classification, 
                     DeCAF_graded = guess)
  rownames(final) = make.names(colnames(data), unique = any(table(colnames(data)) > 1) )
  
  return(final)
}
