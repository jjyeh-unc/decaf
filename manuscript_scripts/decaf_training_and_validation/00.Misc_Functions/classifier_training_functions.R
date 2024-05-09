### A set of functions used to create and evaluate classifier

fitfunc = function(classifier, train_sub, class, skip = F,  alpha = 1, fold = ncol(train_sub), raneff = F, study = NULL,npen = 50, penalty = "MCP", settonull = F){
  indmat = matrix(-1, ncol(train_sub), nrow(classifier$TSPs))
  for(i in 1:nrow(classifier$TSPs)){
    p1 = which(rownames(train_sub) == classifier$TSPs[i,1])
    p2 = which(rownames(train_sub) == classifier$TSPs[i,2])
    indmat[,i] = (train_sub[p1,] > train_sub[p2,])^2
  }
  colnames(indmat) = paste("indmat", 1:ncol(indmat), sep = "")
  data = data.frame(indmat, class = factor(class))
  form = "class"
  add = rep("+", ncol(indmat))
  add[1] = "~"
  for(i in 1:ncol(indmat)) form = sprintf("%s %s %s", form, add[i], colnames(indmat)[i])
  form = as.formula(form)
  
  if(skip == T){
    if(raneff == F){
      library(ncvreg)
      class = as.numeric(class)-1
      cv =  cv.ncvreg(y = class, X =indmat, family = "binomial", penalty = penalty, nfolds = fold , alpha = alpha)
      fit = ncvreg(y = class, X =indmat, family = "binomial", lambda = cv$lambda.min, penalty = penalty, alpha = alpha)
      fit$pe = cv$pe[cv$min]
    }else{
      library(glmmLasso)
      if(settonull == T){
        rnd = NULL
      }else{
        rnd = list(study = ~ 1)
      }
      class = as.numeric(class)-1
      data2  = data.frame( scale(indmat,center = T, scale = T), class = class, study = study)
      classifier$data2 = data2
      pen = exp(seq(log(10), log(1000),,npen))
      bic = coef = rep(10^10, length(pen))
      for(i in 1:length(pen)){
        out = glmmLasso(fix= form , rnd=rnd, data = data2, lambda = pen[i], family = binomial(),final.re=TRUE)
        bic[i] = out$bic
        coef[i] = sum(abs(out$coef[-1]) > 0)
        print(i)
      }
      fit = glmmLasso(fix= form , rnd=rnd, data = data2, lambda = pen[which.min(bic)], family = binomial(),final.re=TRUE)
      
    }
  }else{
    
    fit = glm(form, data = data, family = binomial())
  }
  
  classifier$fit = fit
  
   return(classifier)
}





create.classif = function(dat, classifier, dec= NULL, drop = F, labels = NULL, fit){
  if(!is.null(dec)){
     }else{
    tp <- SWAP.KTSP.Classify(inputMat = dat, classifier = classifier)
  }
  
  if(is.null(dec)){
    p2 = as.numeric(as.matrix(tp)) - 1   
  }else{
    if(any(names(classifier) == "skip")){
      indmat = matrix(-1, ncol(dat), nrow(dat))
      for(i in 1:ncol(indmat)){
        indmat[,i] = (dat[i,] > median(dat[i,]))^2
      }
      colnames(indmat) = rownames(dat)
      data = data.frame(indmat)
      p2 = predict(fit, data)$pred.prob[,2]
    }else{
      indmat = matrix(-1, ncol(dat), nrow(classifier$TSPs))
      for(i in 1:nrow(classifier$TSPs)){
        p1 = which(rownames(dat) == classifier$TSPs[i,1])
        p2 = which(rownames(dat) == classifier$TSPs[i,2])
        if(length(p1) == 0 | length(p2) == 0){
          indmat[i,] = 0
          print(paste("skipping",classifier$TSPs[i,]))
        }else{
          indmat[,i] = (dat[p1,] > dat[p2,])^2
        }
      }
      colnames(indmat) = paste("indmat", 1:ncol(indmat), sep = "")
      
      if(any(names(fit) == "a0")){
        X=cbind(rep(1, nrow(indmat)), indmat)
        p2 = exp(X%*%c(fit$a0,as.numeric(fit$beta)))/(1+exp(X%*%c(fit$a0,as.numeric(fit$beta))))
        data = indmat
      }else if(any(names(fit) == "convex.min")){
        X=cbind(rep(1, nrow(indmat)), indmat)
        p2 = exp(X%*%c(fit$beta))/(1+exp(X%*%c(fit$beta)))
        data = indmat
      }else if(any(names(fit) == "Deltamatrix")){
        X=cbind(rep(1, nrow(indmat)), indmat)
        p2 = exp(X%*%c(fit$coef))/(1+exp(X%*%c(fit$coef)))
        data = indmat
      }else{
        data = data.frame(indmat)
        p2 = predict(fit, data, type = "response")
      }
    }
    tp = (p2 > .5)^2
    names(tp) = names(p2) = colnames(dat)
    if(ncol(dat) != nrow(data)){
      print(dim(dat))
      print(dim(data))
      print(dim(indmat))
      print(dim(train_sub3))
    }
  }
  t = list(class = tp, predprob = p2, dat = dat,classifier = classifier, dec = dec, labels = labels)
  return(t)
}



dec = function(x){
  load("temp")
  data = data.frame(matrix((x)^2, nrow = 1))
  colnames(data) = paste("indmat", 1:(ncol(fit$data)-1), sep = "")
  p = predict.glm(fit, newdata = data, type = "response")
  return(p > .5)
}

