#' Support functions
#'
#' This function define the different techniques for Survival, Binary or Continuous
#'
penaltytech = function(df, outvar, modeltype = "regression",...){

  set.seed(1)
  if(modeltype=="survival"){Y=cbind(time=df[,outvar[1]], status=df[,outvar[2]]); family = "cox"}
  else if(modeltype=="logistic"){Y=as.matrix(df[,outvar]); family = "binomial"}
  else{Y=as.matrix(df[,outvar]); family = "gaussian"}

  if(ncol(df)-length(outvar) > 1){
    inputvarnames = names(df)[!names(df) %in% outvar]
    inputclass = sapply(inputvarnames, function(x) {class(df[,x])})
    if(any(inputclass == "factor")){dfmat = model.matrix( ~ ., df[,!names(df) %in% outvar])[, -1]}
    else{dfmat = data.matrix(df[,!names(df) %in% outvar])}
    }
  else{
    tempdf = data.frame(df[,!names(df) %in% outvar])
    names(tempdf) = names(df)[!names(df) %in% outvar]
    dfmat = model.matrix( ~ ., tempdf)
    if(ncol(dfmat)>2){dfmat = dfmat[,-1]}
  }
  arglist= list(...)

  if(any(arglist$penalty.factor %in% c("adaptive")) | any(!arglist$alpha %in% c(0,1))){
    if(any(arglist$penalty.factor %in% c("adaptive"))){
      if(length(arglist$penalty.factor) == 1){arglist$penalty.factor = rep(1,ncol(df)-length(outvar))}
      else{arglist$penalty.factor = arglist$penalty.factor[-1]}

      #Run Ridge
      cv.ridge <- glmnet::cv.glmnet(dfmat, Y, alpha=0, standardize=F, family=family, penalty.factor = arglist$penalty.factor)
      if(modeltype == "survival"){
        ## Using gamma = 1
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[1:(ncol(dfmat)), 1]))^0.25
      }
      else{
        ## Using gamma = 1
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(dfmat)+1), 1]))^0.25
      }
      dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
      arglist$penalty.factor[which(arglist$penalty.factor != 0)]=dif_wt[which(arglist$penalty.factor != 0)]
    }
    if(any(!arglist$alpha %in% c(0,1))){
      if(length(arglist$penalty.factor) == 0){arglist$penalty.factor = rep(1,ncol(df)-length(outvar))}
      alpha_val = sapply(1:9, function(x){
                              set.seed(1)
                              fit.elnet.cv <- glmnet::cv.glmnet(dfmat, Y, alpha=x/10, standardize=F, family=family, penalty.factor = arglist$penalty.factor) # get optimum lambda
                              cvmse = fit.elnet.cv$cvm[which(fit.elnet.cv$lambda ==fit.elnet.cv$lambda.min)][1]
                              cvmse
      })
      sel_alpha = which(alpha_val == min(alpha_val, na.rm = T))/10
      arglist$alpha = sel_alpha
    }

    # Perform glmnet
    set.seed(1)
    cvfit = glmnet::cv.glmnet(dfmat, Y, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor) # get optimum lambda

    lambda.1se=cvfit$lambda.1se
    lambda.min=cvfit$lambda.min
    nonzerocoef = which(coef(cvfit, s=lambda.1se)[,1] != 0)
    if(length(nonzerocoef) <= 1){
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
    }
    else{
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
    }
  }
  else{

    # Perform glmnet
    fit = glmnet::cv.glmnet(dfmat, Y, standardize=F, family= family, ...) # get optimum lambda
    lambda.1se=fit$lambda.1se
    lambda.min=fit$lambda.min
    nonzerocoef = which(coef(fit, s=lambda.1se)[,1] != 0)

    if(length(nonzerocoef) <= 1){
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, standardize=F, family= family, ...)
    }
    else{
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, standardize=F, family= family, ...)
    }
  }

  fit = model
  return(fit)
}

penaltytech_int = function(df, outvar, modeltype = "regression",...){

  set.seed(1)
  if(modeltype=="survival"){Y=cbind(time=df[,outvar[1]], status=df[,outvar[2]]); family = "cox"}
  else if(modeltype=="logistic"){Y=as.matrix(df[,outvar]); family = "binomial"}
  else{Y=as.matrix(df[,outvar]); family = "gaussian"}

  if(ncol(df)-length(outvar) > 1){
    inputvarnames = names(df)[!names(df) %in% outvar]
    inputclass = sapply(inputvarnames, function(x) {class(df[,x])})
    dfmat = model.matrix( ~ .*., df[,!names(df) %in% outvar])[, -1]
    # if(any(inputclass == "factor")){dfmat = model.matrix( ~ .*., df[,!names(df) %in% outvar])[, -1]}
    # else{dfmat = data.matrix(df[,!names(df) %in% outvar])}
  }
  else{
    tempdf = data.frame(df[,!names(df) %in% outvar])
    names(tempdf) = names(df)[!names(df) %in% outvar]
    dfmat = model.matrix( ~ ., tempdf)
    if(ncol(dfmat)>2){dfmat = dfmat[,-1]}
  }
  arglist= list(...)

  if(any(arglist$penalty.factor %in% c("adaptive")) | any(!arglist$alpha %in% c(0,1))){
    if(any(arglist$penalty.factor %in% c("adaptive"))){
      #Run Ridge
      cv.ridge <- glmnet::cv.glmnet(dfmat, Y, alpha=0, standardize=F, family=family)
      if(modeltype == "survival"){
        ## Using gamma = 1
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[1:(ncol(dfmat)), 1]))^0.25
      }
      else{
        ## Using gamma = 1
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(dfmat)+1), 1]))^0.25
      }
      dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
      arglist$penalty.factor=dif_wt
    }
    if(any(!arglist$alpha %in% c(0,1))){
      if(length(arglist$penalty.factor) == 0){arglist$penalty.factor = rep(1,ncol(df)-length(outvar))}
      alpha_val = sapply(1:9, function(x){
        set.seed(1)
        fit.elnet.cv <- glmnet::cv.glmnet(dfmat, Y, alpha=x/10, standardize=F, family=family, penalty.factor = arglist$penalty.factor) # get optimum lambda
        cvmse = fit.elnet.cv$cvm[which(fit.elnet.cv$lambda ==fit.elnet.cv$lambda.min)][1]
        cvmse
      })
      sel_alpha = which(alpha_val == min(alpha_val, na.rm = T))/10
      arglist$alpha = sel_alpha
    }

    # Perform glmnet
    set.seed(1)
    cvfit = glmnet::cv.glmnet(dfmat, Y, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor) # get optimum lambda

    lambda.1se=cvfit$lambda.1se
    lambda.min=cvfit$lambda.min
    nonzerocoef = which(coef(cvfit, s=lambda.1se)[,1] != 0)
    if(length(nonzerocoef) <= 1){
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
    }
    else{
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
    }
  }
  else{

    # Perform glmnet
    fit = glmnet::cv.glmnet(dfmat, Y, standardize=F, family= family, ...) # get optimum lambda
    lambda.1se=fit$lambda.1se
    lambda.min=fit$lambda.min
    nonzerocoef = which(coef(fit, s=lambda.1se)[,1] != 0)

    if(length(nonzerocoef) <= 1){
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, standardize=F, family= family, ...)
    }
    else{
      model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, standardize=F, family= family, ...)
    }
  }

  fit = model
  return(fit)
}

#'@export
lasso = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){

  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]

  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 1)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1)}
  }

  return(fit)
}
#'@export
alasso = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){
  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]
  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    penalty.factor = c("adaptive", penalty.factor)

    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = "adaptive")}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = "adaptive")}
  }

  return(fit)
}

#'@export
ridge = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){
  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]
  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0)}
  }

  return(fit)
}
#'@export
aridge = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){
  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]
  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    penalty.factor = c("adaptive", penalty.factor)

    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = "adaptive")}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = "adaptive")}
  }

  return(fit)
}

#'@export
enet = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){
  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]
  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1)}
  }

  return(fit)
}
#'@export
aenet = function(df, outvar, modeltype = "regression", interaction = F, covs = NA){
  # df = df[,apply(df, 2, var, na.rm=TRUE) != 0]
  if(!is.na(covs)){
    covindex = which(names(df) %in% covs)
    penalty.factor = rep(1,ncol(df) - length(outvar))
    penalty.factor[covindex] = 0
    penalty.factor = c("adaptive", penalty.factor)

    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = penalty.factor)}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = penalty.factor)}
  }
  else{
    if(interaction){fit = penaltytech_int(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = "adaptive")}
    else{fit = penaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = "adaptive")}
  }
  return(fit)
}

fastpenaltytech = function(df, outvar, modeltype = "regression",...){
  set.seed(1)
  if(modeltype=="survival"){Y=cbind(time=df[,outvar[1]], status=df[,outvar[2]]); family = "cox"}
  else if(modeltype=="logistic"){Y=as.matrix(df[,outvar]); family = "binomial"}
  else{Y=as.matrix(df[,outvar]); family = "gaussian"}

  if((ncol(df)-length(outvar)) > 1){
    inputvarnames = names(df)[!names(df) %in% outvar]
    inputclass = sapply(inputvarnames, function(x) {class(df[,x])})
    if(any(inputclass == "factor")){dfmat = model.matrix( ~ ., df[,!names(df) %in% outvar])[, -1]}
    else{dfmat = data.matrix(df[,!names(df) %in% outvar])}
  }
  else{
    dfmat = model.matrix( ~ ., df[,!names(df) %in% outvar])}
  arglist= list(...)

  cond1 = any(names(arglist) %in% "cv")
  cond2 = any(names(arglist) %in% "hyperpara")
  cond3 = arglist$hyperpara == T

  if(cond1){cond = T}
  else if(cond2){if(cond3){cond = T}else{cond = F}}
  else{cond = F}

  if(cond){
    cv = ifelse(is.null(arglist$cv), 10,arglist$cv)
    if(any(arglist$penalty.factor %in% c("adaptive")) | any(!arglist$alpha %in% c(0,1))){

      if(any(arglist$penalty.factor %in% c("adaptive"))){
        #Run Ridge
        cv.ridge <- glmnet::cv.glmnet(dfmat, Y, alpha=0, standardize=F, family=family, nfolds = cv)
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(dfmat)+1), 1]))^0.25 ## Using gamma = 1
        dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
        arglist$penalty.factor=dif_wt
      }
      if(any(!arglist$alpha %in% c(0,1))){
        if(length(arglist$penalty.factor) == 0){arglist$penalty.factor = rep(1,ncol(df)-length(outvar))}
        alpha_val = sapply(1:9, function(x){
          set.seed(1)
          fit.elnet.cv <- glmnet::cv.glmnet(dfmat, Y, alpha=x/10, standardize=F, family=family, penalty.factor = arglist$penalty.factor, nfolds = cv) # get optimum lambda
          cvmse = fit.elnet.cv$cvm[which(fit.elnet.cv$lambda ==fit.elnet.cv$lambda.min)][1]
          cvmse
        })
        sel_alpha = which(alpha_val == min(alpha_val, na.rm = T))/10
        arglist$alpha = sel_alpha
      }

      # Perform glmnet
      set.seed(1)
      cvfit = glmnet::cv.glmnet(dfmat, Y, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor, nfolds = cv) # get optimum lambda

      lambda.1se=cvfit$lambda.1se
      lambda.min=cvfit$lambda.min
      nonzerocoef = which(coef(cvfit, s=lambda.1se)[,1] != 0)
      if(length(nonzerocoef) <= 1){
        model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
      }
      else{
        model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, alpha=arglist$alpha, standardize=F, family= family, penalty.factor = arglist$penalty.factor)
      }
    }
    else{
      argvec = c(...)
      argvec = argvec[which(names(argvec) != "cv")]
      # Perform glmnet
      fit = do.call(glmnet::cv.glmnet,c(argvec, list(x = dfmat, y = Y, standardize=F, family= family, nfolds = cv)))
      # fit = glmnet::cv.glmnet(dfmat, Y, standardize=F, family= family, nfolds = cv, ...) # get optimum lambda

      lambda.1se=fit$lambda.1se
      lambda.min=fit$lambda.min
      nonzerocoef = which(coef(fit, s=lambda.1se)[,1] != 0)

      if(length(nonzerocoef) <= 1){
        model = do.call(glmnet::glmnet,c(argvec, list(x = dfmat, y = Y, lambda = lambda.min, standardize=F, family= family)))
        # model = glmnet::glmnet(dfmat, Y, lambda = lambda.min, standardize=F, family= family, ...)
      }
      else{
        model = do.call(glmnet::glmnet,c(argvec, list(x = dfmat, y = Y, lambda = lambda.1se, standardize=F, family= family)))
        # model = glmnet::glmnet(dfmat, Y, lambda = lambda.1se, standardize=F, family= family, ...)
      }
    }
  }
  else{
    if(any(arglist$penalty.factor %in% c("adaptive")) | any(!arglist$alpha %in% c(0,1))){

      if(any(arglist$penalty.factor %in% c("adaptive"))){
        #Run Ridge
        cv.ridge <- glmnet::cv.glmnet(dfmat, Y, alpha=0, standardize=F, family=family)
        dif_wt <- 1/abs(matrix(coef(cv.ridge, s=cv.ridge$lambda.min)[2:(ncol(dfmat)+1), 1]))^0.25 ## Using gamma = 1
        dif_wt[dif_wt[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
        arglist$penalty.factor=dif_wt
      }

      # Perform glmnet
      set.seed(1)
      argvec = unlist(arglist)
      argvec = argvec[which(names(argvec) != "cv")]

      # Perform glmnet
      model = do.call(glmnet::glmnet,c(argvec, list(x = dfmat, y = Y, standardize=F, family= family)))
    }
    else{
      argvec = c(...)
      argvec = argvec[which(names(argvec) != "cv")]
      # Perform glmnet
      model = do.call(glmnet::glmnet,c(argvec, list(x = dfmat, y = Y, standardize=F, family= family)))
    }
  }

  fit = model
  return(fit)
}

#'@export
flasso = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, ...)
  return(fit)
}
#'@export
falasso = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 1, penalty.factor = "adaptive",...)
  return(fit)
}
#'@export
fridge = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0,...)
  return(fit)
}
#'@export
faridge = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0, penalty.factor = "adaptive",...)
  return(fit)
}

#'@export
fenet = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1,...)
  return(fit)
}
#'@export
faenet = function(df, outvar, modeltype = "regression",...){
  fit = fastpenaltytech(df = df, outvar = outvar, modeltype = modeltype, alpha = 0.1, penalty.factor = "adaptive",...)
  return(fit)
}
