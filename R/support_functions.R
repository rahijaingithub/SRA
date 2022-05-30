#' Support functions
#'
#' This function run the different techniques for Survival, Binary or Continuous
modelfunction = function(seed = 1, traindf, outvar, modeltype, technique, hyperpara = NA){
  set.seed(seed)
  if(technique == "rf"){
    fit = rftech(df = traindf, outvar, modeltype)

    # if(modeltype != "survival"){
    #   f = as.formula(paste(outvar, "~."))
    #   fit <- randomForest::randomForest(f, data=traindf)
    # }
    # else{
    #   f = as.formula(paste("Surv(SurvTime, SurvEvent)", "~."))
    #   fit <- randomForestSRC::rfsrc(f, traindf, ntree = 500, samptype = "swr")
    # }
  }
  else if(technique %in% c("lasso", "LASSO", "ridge", "enet", "alasso", "aenet", "aridge")){
    modellist = c(lasso = lasso, LASSO = lasso, ridge = ridge, enet = enet, alasso = alasso, aenet = aenet, aridge = aridge)
    fit = modellist[[technique]](df = traindf, outvar, modeltype)
  }
  else{fit = NA}
  return(fit)
}

#' This function evaluates the performance of different techniques
accuracy = function(pred, actual){
  acc = MLmetrics::Accuracy(pred, actual)
  c_accu = 1-acc
  perf = max(c(acc, 1-acc), na.rm = T)
  return(perf)
}
prauc = function(pred, actual){ # It treats one as event/postive/target feature and 0 as noise/noevent/negative
  perf = PRROC::pr.curve(scores.class0 = pred[which(actual ==1)], scores.class1 = pred[which(actual ==0)])$auc.integral
  return(perf)
}
#'@export
perffunct = function(modeltype, outvar, testdf, pred = NA, technique = "rf", model){

  # Check if prediction needs to be done or not
  if(is.na(pred)){
    if(technique %in% c("rf", "autorf")){
      if(modeltype == "survival"){
        predvar = predict(model, newdata = testdf)
        pred = predvar$predicted
      }
      else{pred = predict(model, newdata = testdf)}
      }
    else if(technique %in% c("ann", "autoann")){
      Predict= neuralnet::compute(model,testdf)
      if(modeltype == "regression"){pred <- Predict$net.result}
      else if(modeltype == "logistic"){
        prob <- Predict$net.result[,2]
        pred <- ifelse(prob>0.5, 1, 0)
        }
      else{pred = pred <- Predict$net.result}
    }
    else{
      if(modeltype == "survival"){f = as.formula(paste("Surv(SurvTime, SurvEvent)~."))}
      else{f = as.formula(paste(outvar,"~."))}

      # inputvarnames = names(testdf)[!names(testdf) %in% outvar]
      # inputclass = sapply(inputvarnames, function(x) {class(testdf[,x])})
      # if(any(inputclass == "factor")){testmat = model.matrix( ~ ., testdf[,!names(testdf) %in% outvar])[, -1]}
      #
      # else{testmat = data.matrix(testdf[,!names(testdf) %in% outvar])}

      if(ncol(testdf)-length(outvar) > 1){
        inputvarnames = names(testdf)[!names(testdf) %in% outvar]
        inputclass = sapply(inputvarnames, function(x) {class(testdf[,x])})
        if(any(inputclass == "factor")){testmat = model.matrix( ~ ., testdf[,!names(testdf) %in% outvar])[, -1]}
        else{testmat = data.matrix(testdf[,!names(testdf) %in% outvar])}
      }
      else{
        tempdf = data.frame(testdf[,!names(testdf) %in% outvar])
        names(tempdf) = names(testdf)[!names(testdf) %in% outvar]
        testmat = model.matrix( ~ ., tempdf)
        if(ncol(testmat)>2){testmat = testmat[,-1]}
      }

      # # testmat = model.matrix(f, testdf)[,-1]
      # testmat = data.matrix(testdf[,!names(testdf) %in% outvar])
      if(modeltype == "logistic"){pred = predict(model, newx = testmat, type = "class")}
      else{pred = predict(model, newx = testmat)}
      if(class(pred) %in% c("numeric", "integer")){pred = pred}
      else{pred = pred[,1]}
    }
  }

  # Get the performance
  if(modeltype == "survival"){
    cind = glmnet::Cindex(pred = pred, y = cbind(time = testdf[,"SurvTime"], status = testdf[,"SurvEvent"]))
    c_cind = 1-cind
    perf = max(c(cind, c_cind), na.rm = T)
    # if(perf < 0){perf = NA}
  }
  else if(modeltype == "regression"){
    perf = MLmetrics::RMSE(pred, testdf[,outvar])
    perf = 1/perf
  }
  else{
    perf = accuracy(pred = pred, actual = testdf[,outvar])
  }

  return(perf)
}
#'@export
perffunct_int = function(modeltype, outvar, testdf, pred = NA, technique = "rf", model){

  # Check if prediction needs to be done or not
  if(is.na(pred)){
    if(technique %in% c("rf", "autorf")){
      if(modeltype == "survival"){
        predvar = predict(model, newdata = testdf)
        pred = predvar$predicted
      }
      else{pred = predict(model, newdata = testdf)}
    }
    else if(technique %in% c("ann", "autoann")){
      Predict= neuralnet::compute(model,testdf)
      if(modeltype == "regression"){pred <- Predict$net.result}
      else if(modeltype == "logistic"){
        prob <- Predict$net.result[,2]
        pred <- ifelse(prob>0.5, 1, 0)
      }
      else{pred = pred <- Predict$net.result}
    }
    else{
      if(modeltype == "survival"){f = as.formula(paste("Surv(SurvTime, SurvEvent)~."))}
      else{f = as.formula(paste(outvar,"~."))}

      # inputvarnames = names(testdf)[!names(testdf) %in% outvar]
      # inputclass = sapply(inputvarnames, function(x) {class(testdf[,x])})
      # if(any(inputclass == "factor")){testmat = model.matrix( ~ ., testdf[,!names(testdf) %in% outvar])[, -1]}
      #
      # else{testmat = data.matrix(testdf[,!names(testdf) %in% outvar])}

      if(ncol(testdf)-length(outvar) > 1){
        inputvarnames = names(testdf)[!names(testdf) %in% outvar]
        inputclass = sapply(inputvarnames, function(x) {class(testdf[,x])})
        testmat = model.matrix( ~ .*., testdf[,!names(testdf) %in% outvar])[, -1]
        # if(any(inputclass == "factor")){testmat = model.matrix( ~ ., testdf[,!names(testdf) %in% outvar])[, -1]}
        # else{testmat = data.matrix(testdf[,!names(testdf) %in% outvar])}
      }
      else{
        tempdf = data.frame(testdf[,!names(testdf) %in% outvar])
        names(tempdf) = names(testdf)[!names(testdf) %in% outvar]
        testmat = model.matrix( ~ ., tempdf)
        if(ncol(testmat)>2){testmat = testmat[,-1]}
      }

      # # testmat = model.matrix(f, testdf)[,-1]
      # testmat = data.matrix(testdf[,!names(testdf) %in% outvar])
      if(modeltype == "logistic"){pred = predict(model, newx = testmat, type = "class")}
      else{pred = predict(model, newx = testmat)}
      if(class(pred) %in% c("numeric", "integer")){pred = pred}
      else{pred = pred[,1]}
    }
  }

  # Get the performance
  if(modeltype == "survival"){
    cind = glmnet::Cindex(pred = pred, y = cbind(time = testdf[,"SurvTime"], status = testdf[,"SurvEvent"]))
    c_cind = 1-cind
    perf = max(c(cind, c_cind), na.rm = T)
    # if(perf < 0){perf = NA}
  }
  else if(modeltype == "regression"){
    perf = MLmetrics::RMSE(pred, testdf[,outvar])
    perf = 1/perf
  }
  else{
    perf = accuracy(pred = pred, actual = testdf[,outvar])
  }

  return(perf)
}

#' This function evaluates the different models used in GA function
modelanalysis = function(df_train, df_test, outvar, modeltype, technique, hyperpara = NA){

  # Run the model
  fit = modelfunction(seed = 1, traindf=df_train, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = hyperpara)
  # Predict the Test Data
  pred = predict(fit, newdata = df_test)
  # Model Performance
  if(modeltype == "survival"){perf = perffunct(modeltype, outvar, testdf = df_test, pred$predicted)}
  else{perf = perffunct(modeltype, outvar, testdf=df_test, pred)}
  return(perf)
}

penal_modelanalysis = function(df_train, df_test, outvar, modeltype, technique, hyperpara = NA, totalvar = 200){

  # Run the model
  fit = modelfunction(seed = 1, traindf=df_train, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = hyperpara)

  # Predict the Test Data
  pred = predict(fit, newdata = df_test)

  # Model Performance
  if(modeltype == "survival"){
    perf = perffunct(modeltype, outvar, testdf = df_test[complete.cases(df_test),], pred$predicted);
    feat = ncol(df_train)-2}
  else{perf = perffunct(modeltype, outvar, testdf=df_test, pred); feat = ncol(df_train)-1}
  if(!is.na(perf)){perf = perf - (1*feat/totalvar)}

  return(perf)
}

#' Perform GA
gamodel = function(f, gatype = c("binary", "real"), maxfeat = NA, traindf, outvar,  ...){

  # Get the index for outcome feature
  outindex_train = which(names(traindf) %in% outvar)

  # Define GA parameters
  if(is.na(maxfeat)){maxfeat = min(200, ncol(traindf)-length(outindex_train))}
  if(gatype == "binary"){
    nBits = ncol(traindf)-length(outindex_train);
    type = "binary";
    lower = NA; upper = NA}
  else{
    lower = rep(1,maxfeat);
    upper = rep(ncol(traindf)-length(outindex_train), maxfeat);
    type = "real-valued";
    nBits = NA}

  arg = list(...)
  defaultgapara = list(popSize = 10, pmutation = 0.3, pcrossover = 0.8, maxiter = 100, run = 1, seed =3)
  missname = setdiff(names(defaultgapara), names(arg))
  arg[missname] = lapply(defaultgapara[missname],c)

  # Run GA
  GA = GA::ga(type = type, fitness = f, nBits = nBits, lower = lower, upper = upper, ...)
  return(GA)
}

#' Perform CV
cvfun_penalty = function(traindf, kfold, outvar, colnum, modeltype, technique, totalvar){
  outindex_train = which(names(traindf) %in% outvar)

  set.seed(1)
  #Randomly shuffle the data
  traindf<-traindf[sample(nrow(traindf)),]

  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(traindf)),breaks=kfold,labels=FALSE)

  #Perform 10 fold cross validation
  cvfs = sapply(1:kfold, function(i) {
    sampleid = which(folds == i)
    cv_traindf = traindf[-sampleid,c(colnum, outindex_train)]
    cv_testdf = traindf[sampleid,c(colnum, outindex_train)]

    perf = penal_modelanalysis(df_train = cv_traindf, df_test = cv_testdf, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = NA, totalvar = totalvar)
  })
  return(cvfs)
}

#' Perform boot
bootfun_penalty = function(traindf, kfold, outvar, colnum, modeltype, technique, totalvar){
  outindex_train = which(names(traindf) %in% outvar)

  if(length(outvar) >1){
    trydf = sapply(traindf[,outvar], function(x) {length(unique(x))})
    stratvar = names(trydf)[trydf == 2]}
  else{
    stratvar = outvar
    trydf = length(unique(traindf[,outvar]))
    if(trydf == 2){stratvar}
    else{stratvar = NULL}
  }

  bootfs = sapply(1:kfold, function(i) {
    set.seed(i)
    if(is.null(stratvar)){
      rowindex= sampling::strata(traindf, stratanames = stratvar, size = nrow(traindf), method="srswr")
    }
    else{rowindex= sampling::strata(traindf, stratanames = stratvar, size = table(traindf[,stratvar]), method="srswr")}
    sampleid = rowindex$ID_unit
    bt_traindf = traindf[sampleid,c(colnum, outindex_train)]
    bt_testdf = traindf[-sampleid,c(colnum, outindex_train)]
    perf = penal_modelanalysis(df_train = bt_traindf, df_test = bt_testdf, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = NA, totalvar = totalvar)
  })

  return(bootfs)
}

#' Perform CV
cvfun = function(traindf, kfold, outvar, colnum, modeltype, technique){
  outindex_train = which(names(traindf) %in% outvar)

  set.seed(1)
  #Randomly shuffle the data
  traindf<-traindf[sample(nrow(traindf)),]

  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(traindf)),breaks=kfold,labels=FALSE)

  #Perform 10 fold cross validation
  cvfs = sapply(1:kfold, function(i) {
    sampleid = which(folds == i)
    cv_traindf = traindf[-sampleid,c(colnum, outindex_train)]
    cv_testdf = traindf[sampleid,c(colnum, outindex_train)]

    perf = penal_modelanalysis(df_train = cv_traindf, df_test = cv_testdf, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = NA)
  })
  return(cvfs)
}

#' Perform boot
bootfun = function(traindf, kfold, outvar, colnum, modeltype, technique){
  outindex_train = which(names(traindf) %in% outvar)
  if(length(outvar) >1){
    trydf = sapply(traindf[,outvar], function(x) {length(unique(x))})
    stratvar = names(trydf)[trydf == 2]}
  else{
    stratvar = outvar
    trydf = length(unique(traindf[,outvar]))
    if(trydf == 2){stratvar}
    else{stratvar = NULL}
    }

  bootfs = sapply(1:kfold, function(i) {

    set.seed(i)

    if(is.null(stratvar)){
      rowindex= sampling::strata(traindf, stratanames = stratvar, size = nrow(traindf), method="srswr")
    }
    else{rowindex= sampling::strata(traindf, stratanames = stratvar, size = table(traindf[,stratvar]), method="srswr")}

    sampleid = rowindex$ID_unit
    bt_traindf = traindf[sampleid,c(colnum, outindex_train)]
    bt_testdf = traindf[-sampleid,c(colnum, outindex_train)]
    perf = modelanalysis(df_train = bt_traindf, df_test = bt_testdf, outvar = outvar, modeltype = modeltype, technique = technique, hyperpara = NA)
  })

  return(bootfs)
}

#'Normalize data
normalized = function(x){y = (x-min(x))/(max(x)-min(x)); return(y)}

#'Output data process
loopdf = function(df, tech){
  # This function extracts feature importance and normalize it between 0 and 1 with 1 being maximum importance.

  if(any(tech %in% c("rf", "rfauto"))){
    if(is.null(df$importance)){
      varimp = randomForestSRC::vimp(df)$importance
      varimp_norm = exp(varimp) # take exponential to convert negative to positive.
      selfeat = data.frame(variable = names(varimp_norm), importance=varimp_norm/max(varimp_norm), stringsAsFactors = F)
    }
    else{
      selfeat = data.frame(variable = rownames(df$importance), importance=df$importance/max(df$importance), stringsAsFactors = F)
    }

  }
  else if(any(tech %in% c("ann", "autoann"))){
    resimp = NeuralNetTools::olden(df, bar_plot=F)
    selfeat = data.frame(variable = rownames(resimp), importance=resimp$importance/max(abs(resimp$importance)), stringsAsFactors = F)
  }
  else{ # For LASSO based methods
    if(any(class(df) %in% "glmnet")){df = glmnet::coef.glmnet(df)}
    selfeat = data.frame(variable = rownames(df)[-1], importance=df[-1,1]/max(abs(df[-1,1])), stringsAsFactors = F)
  }
  names(selfeat) = c("variable", "importance")
  selfeat$importance[is.na(selfeat$importance) | is.nan(selfeat$importance)] = 0
  return(selfeat)
}
outdat = function(tech, resdf, pretech, martarget = 10, algo = c("AIFS", "Std"), perfmetric = c("AUC", "Accuracy", "F1")){

  resdf[,2] = abs(resdf[,2])
  resdf[is.nan(resdf[,2]) | is.na(resdf[,2]),2] = 0

  if(any(tech %in% c("rf", "ridge", "aridge"))){
    resdf = resdf[order(resdf[,2], decreasing = T),]
    resdf$predlab = 0
    resdf$predlab[1:martarget] = 1
  }
  else if(any(tech %in% c("km"))){
    resdf = resdf[order(resdf[,2], decreasing = T),]
    resdf$predlab = 0
    resdf$predlab[resdf[,3] == min(resdf[,3]) & resdf[,2] != 0] = 1
  }
  else{
    resdf[,2] = abs(resdf[,2])
    resdf = resdf[order(resdf[,2], decreasing = T),]
    resdf$predlab = 0
    resdf$predlab[resdf[,2] !=0] = 1
  }

  tarfeat = paste("X",1:martarget, sep="")
  resdf$truelab = 0
  resdf$truelab[resdf$variable %in% tarfeat] = 1

  if(perfmetric == "AUC"){resdf$perf = MLmetrics::AUC(y_pred = resdf$predlab, y_true =  as.factor(resdf$truelab))}
  else if(perfmetric == "Accuracy"){resdf$perf = MLmetrics::Accuracy(y_pred = resdf$predlab, y_true = as.factor(resdf$truelab))}
  else{
    if(all(resdf$predlab != 1)){resdf$perf = 0}
    else if(all(resdf$predlab != 0)){
      precision = MLmetrics::Precision(y_pred = resdf$predlab, y_true = as.factor(resdf$truelab), positive = 1)
      recall = 1
      resdf$perf = 2*precision*recall/(precision + recall)}
    else{resdf$perf = MLmetrics::F1_Score(y_pred = resdf$predlab, y_true = as.factor(resdf$truelab), positive = 1)}
    resdf$perf[is.nan(resdf$perf)] = 0
    }
  resdf$type = paste(algo,pretech,tech, sep = "_")

  return(resdf)
}
colcint = function(column){
  ci = unlist(list(attributes(suppressMessages(miscset:::confint.numeric(column, na.rm = T, level = 0.95)))))
  res = c(lci = ci[2], uci = ci[3])
  return(res)
}


#'Ensemble Sampling
#'@export
bootsample = function(df = traindf, boots = 10, samplesize = 1.0){
  if(is.na(boots)){boots = 100}
  samlist = lapply(1:boots, function(x){ set.seed(x)
                                         sam = sample(rownames(df), nrow(df), replace = T)
                                         sam
                                        })
  return(samlist)
}

#'@export
randomsample = function(df = traindf, boots = 10, samplesize = 0.8){
  if(is.na(boots)){boots = 100}
  samlist = lapply(1:boots, function(x){ set.seed(x)
    sam = sample(rownames(df), floor(nrow(df)*samplesize), replace = F)
    sam
  })
  return(samlist)
}

#'@export
cvsample = function(df = traindf, boots = 10, samplesize = 0.9){
  if(is.na(boots)){boots = 100}

  kfold = floor(1/min((1-samplesize), samplesize))
  shufflenum = floor(boots/kfold)

  samlist = lapply(1:shufflenum, function(x){
                                              set.seed(x)
                                              rowvec = sample(rownames(df), nrow(df), replace = F)
                                              folds = cut(1:nrow(df),breaks=kfold,labels=FALSE)
                                              samfold = lapply(1:kfold, function(y){yfold = rowvec[which(folds!=y)]})
                                              samfold
                                             })

  samlist = purrr::flatten(samlist)
  return(samlist)
}

#'@export
allsample = function(df = traindf, boots = 10, samplesize = 1.0){
  if(is.na(boots)){boots = 100}
  samlist = lapply(1:boots, function(x){ set.seed(x)
    sam = sample(1:nrow(df), floor(nrow(df)*samplesize), replace = F)
    sam
  })
  return(samlist)
}

#'@export
allfeat = function(df = traindf, outvar = outvar, maxfeatsample = ncol(traindf)-length(outvar), boots = 1){
  if(is.na(boots)){boots = 100}
  features = names(df)[!names(df) %in% outvar]
  featlist = lapply(1:boots, function(x){ set.seed(x)
                             sam = sample(features, maxfeatsample, replace = F)
                             sam
                              })
  return(featlist)
}

#'@export
randomfeat = function(df = traindf, outvar = outvar, maxfeatsample = ncol(traindf)-length(outvar), boots = 1){
  if(is.na(boots)){boots = 100}
  features = names(df)[!names(df) %in% outvar]
  featlist = lapply(1:boots, function(x){ set.seed(x)
    sam = sample(features, sample(2:maxfeatsample,1), replace = F)
    sam
  })
  return(featlist)
}

#'Performance Functions
#'@export
predperf = function(traindf, testdf = NA, outvar, modeltype = "regression", martarget =NA, selfeat = NA, pretech, cv_split = 30, savememory = F, Trace = T, parallel = F, smart = F, perfparallel = F, interaction = T){
  if(all(is.na(testdf))){


    if(!is.na(selfeat)){
      if(interaction){
        revfeat = stringi::stri_replace_all_fixed(selfeat, ".", "*")
        f = as.formula(paste("~",paste(revfeat, collapse = "+")))
        dfmat = model.matrix(f, traindf)[,-1]
        dfmat_df = data.frame(dfmat)
        dfmat_df[,outvar] = traindf[,outvar]
        traindf = data.frame(dfmat_df)
      }
      else{traindf = traindf[,c(selfeat, outvar)]}
    }

    # Prepare Model
    fit = ensem_homo(df = traindf, outvar = outvar, modeltype = modeltype, technique = "ridge", sampletype = "cvsample", fstype = "allfeat", maxfeatsample = ncol(traindf)-length(outvar), boots = cv_split, parallel = parallel, perfparallel = perfparallel, Trace=Trace, savememory = savememory, smart=smart)


    # Get Performance for each split
    pp = fit$df

    # Get the feature selection performance
    if(is.na(martarget)){ target = length(selfeat); noise = 0}
    else{
      tarfeat = paste("X",1:martarget, sep="")
      target = length(which(selfeat %in% tarfeat))
      noise = length(selfeat) - target
    }

    # Get overall performance
    cint = colcint(column = pp$perf)
    df = data.frame(pp = mean(pp$perf, na.rm=T), pplci= cint[1], ppuci = cint[2], target =  target, noise = noise, stringsAsFactors = F)
    df$type = paste("ridge",modeltype, pretech, sep = "_")

  }
  else{
    if(length(selfeat) > 0){
      if(!is.na(selfeat)){
        if(interaction){
          revfeat = stringi::stri_replace_all_fixed(selfeat, ".", "*")
          f = as.formula(paste("~",paste(revfeat, collapse = "+")))
          dfmat = model.matrix(f, traindf)[,-1]
          dfmat_df = data.frame(dfmat)
          dfmat_df[,outvar] = traindf[,outvar]
          traindf = data.frame(dfmat_df)

          dfmat = model.matrix(f, testdf)[,-1]
          dfmat_df = data.frame(dfmat)
          dfmat_df[,outvar] = testdf[,outvar]
          testdf = data.frame(dfmat_df)
        }
        else{
          # str(traindf)
          # print(selfeat)
          # print(outvar)
          traindf = traindf[,c(selfeat, outvar)]
          testdf = testdf[,c(selfeat, outvar)]
        }
      }

      # Prepare Model
      fit = emb_mod(df = traindf, outvar = outvar, modeltype = modeltype, technique = "ridge")

      # Get the Prediction Performance
      pp = perffunct(modeltype = modeltype, outvar = outvar, testdf = testdf, technique = "ridge", model = fit)

      # Get the feature selection performance
      if(is.na(martarget)){ target = length(selfeat); noise = 0}
      else{
        tarfeat = paste("X",1:martarget, sep="")
        target = length(which(selfeat %in% tarfeat))
        noise = length(selfeat) - target
      }

      df = data.frame(pp = pp, target =  target, noise = noise, stringsAsFactors = F)
      df$type = paste("ridge",modeltype, pretech, sep = "_")
    }
    else{
      pp = MLmetrics::RMSE(mean(testdf[,outvar]), testdf[,outvar])
      pp = 1/pp
      df = data.frame(pp = pp, target =  0, noise = 0, stringsAsFactors = F)
      df$type = paste("ridge",modeltype, pretech, sep = "_")
    }
  }

  return(df)
}

#'@export
var_organise = function (inputlist, symbol = "_") {
  org_var = lapply(inputlist, function(x) {
    a = stringi::stri_split_fixed(x, symbol)
    ifelse(length(a[[1]]) > 1, paste0(naturalsort::naturalsort(a[[1]]),
                                      collapse = symbol), a[[1]])
  })
  return(unlist(org_var))
}
alltechlist = function(...){
  alltech = expand.grid(..., stringsAsFactors = F)
  return(alltech)
}

modelperf = function(rankdf = finrank, modeltype = "regression", traindf = traindf, testdf = testdf, outvar){
  # Create weighted matrix
  dfmat = lapply(list(traindf,testdf), function(x){
    xmat = x[,!names(x) %in% outvar]
    xmatcol = colnames(xmat)
    rankdf = rankdf[order(match(rankdf$variable,xmatcol)),]
    wt = rankdf$imp/ max(rankdf$imp, na.rm = NA)
    # print(rankdf$imp)
    # print(wt)
    xmat = t(xmat)
    wtmat = t(xmat * wt)
    matsum = rowSums(wtmat)
    df = data.frame(x=matsum, x[,outvar], stringsAsFactors = F)
    names(df) = c("x", outvar)
    if(is.nan(unique(df$x)) | is.na(unique(df$x))){df$x = 0}
    # print(df)
    df
  })
  # str(dfmat)

  # Prepare Model
  if(modeltype == "regression"){
    f = as.formula(paste(outvar,"~x"))
    fit = lm(f, data = dfmat[[1]])
  }
  else if(modeltype == "logistic"){
    f = as.formula(paste(outvar,"~x"))
    fit = glm(f, data = dfmat[[1]], family = "binomial")
  }
  else{
    f = as.formula(paste("Surv(",outvar,")~x"))
    fit = coxph(f, data = dfmat[[1]])
  }

  # Predict Performance
  perf = sapply(dfmat,function(x){
                      perffunct(modeltype = modeltype, outvar = outvar, testdf = x, technique = "reg", model = fit)})
  perf = 1/perf
  names(perf) = c("train", "test")
  return(perf)
}

Correlremoval = function(inputdf, cutoff){
  df = inputdf
  for(i in 1:(ncol(inputdf)-1)){
    # cat(ncol(df), " ")
    x = df[,i]
    y = data.frame(df[,-c(1:i)])
    if(ncol(y) == 0){break}
    dfcor=cor(x=x, y=y, use = "pairwise.complete.obs", method = "spearman")
    dupvar = which(dfcor >= cutoff | dfcor <= -1*cutoff)
    # cat(length(dupvar), " ")
    if(length(dupvar)>0){
      if(ncol(y) == 1){varname = names(df)[i+1]}
      else{varname = names(y)[dupvar]}
      df = df[, !names(df) %in% varname]
    }
    if(ncol(y) == 1){break}
  }
  return(df)
}

