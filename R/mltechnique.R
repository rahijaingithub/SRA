#' Support functions
#'
#' This function define the different techniques for Survival, Binary or Continuous
#'@export
rftech = function(df, outvar, modeltype,...){
  set.seed(1)
  if(modeltype != "survival"){
    f = as.formula(paste(outvar, "~."))
    fit <- randomForest::randomForest(f, data=df, ...)
  }
  else{
    # print(T)
    f = as.formula(paste("Surv(SurvTime, SurvEvent)", "~."))
    fit <- randomForestSRC::rfsrc(f, data = df,...)
    # print("Clear RF")
  }

  return(fit)
}

#'@export
anntech = function(df, outvar, modeltype,...){

  arglist = list(...)
  argvec = c(...)

  # Handle "Hidden" argument
  cond1 = any(grepl("hidden",names(arglist)))
  cond2 = any(names(arglist) %in% "layer")
  if(cond1){if(cond2){cond = T}else{cond = F}}
  else{cond = F}

  if(cond){
    hidden = argvec[grep("hidden",names(argvec))]
    hidden = rep(hidden, arglist$layer);
    hidden = c(as.numeric(hidden[1:arglist$layer]))
  }
  else if(cond1){hidden = argvec[grep("hidden",names(argvec))]}
  else{hidden = 1}

  if(any(names(arglist) %in% "act.fct")){act.fct = arglist$act.fct}else{act.fct = "logistic"}

  argvec = argvec[!(names(argvec) %in% c("layer", "hidden", "act.fct"))]
  argvec = argvec[!grepl("hidden",names(argvec))]
  if(length(argvec) == 0){argvec = NULL}

  numargvec = as.numeric(argvec)
  names(numargvec) = names(argvec)
  argvec = numargvec

  set.seed(1)
  if(modeltype == "regression"){
    n <- names(df)
    f <- as.formula(paste(outvar, "~", paste(n[!n %in% outvar], collapse = " + ")))
    fit <- do.call(neuralnet::neuralnet,c(argvec, list(formula = f, data=df, linear.output = T,
                                                       hidden = hidden, act.fct = act.fct, lifesign = "none")))
  }
  else if(modeltype == "logistic"){
    n <- names(df)
    f <- as.formula(paste(outvar, "~", paste(n[!n %in% outvar], collapse = " + ")))
    fit <- do.call(neuralnet::neuralnet,c(argvec, list(formula = f, data=df, linear.output = F,
                                                       hidden = hidden, act.fct = act.fct)))
  }
  else{
    f = as.formula(paste("Surv(SurvTime, SurvEvent)", "~."))
    fit <- randomForestSRC::rfsrc(f, traindf,...)
  }

  return(fit)
}

#'@export
autorftech = function(df, outvar, modeltype){

  feat = ncol(df)-length(outvar)

  mtry = unique(c(2, max(2, min(floor(feat/3),100)), max(2, min(floor(feat/4),100)), max(2, min(floor(feat/2),100)), max(2, floor(log(feat))),max(2, floor(sqrt(feat))), min(feat,100)))

  ntree = unique(c(250, 500, min(nrow(df),1000), min(ncol(df),1000)))
  if(modeltype == "logistic"){nodesize = c(1:10)}
  else if(modeltype == "survival"){nodesize = c(10:20)}
  else{nodesize = c(5:10)}
  if(modeltype == "survival"){samptype = c("swor", "swr")}
  else{replace = c(T,F)}

  if(modeltype == "survival"){
    para = hyperaparameter(df = df, outvar = outvar, modeltype = modeltype, technique = "rf", mtry = mtry, ntree = ntree, nodesize = nodesize, samptype = samptype)
    # str(para)
    mtry = para$mtry; ntree = para$ntree; nodesize = para$nodesize; samptype = para$samptype
    fit = rftech(df = df, outvar = outvar, modeltype = modeltype, mtry = mtry, ntree = ntree, nodesize = nodesize, samptype = samptype)
  }
  else{
    para = hyperaparameter(df = df, outvar = outvar, modeltype = modeltype, technique = "rf", mtry = mtry, ntree = ntree, nodesize = nodesize, replace = replace)
    mtry = para$mtry; ntree = para$ntree; nodesize = para$nodesize; replace = para$replace
    fit = rftech(df = df, outvar = outvar, modeltype = modeltype, mtry = mtry, ntree = ntree, nodesize = nodesize, replace = replace)
  }
  return(fit)
}

#'@export
autoanntech = function(df, outvar, modeltype){

  feat = ncol(df)-1
  hidden = c(2, 1, 3, max(2, floor(feat/3))) # c(2, max(2, floor(feat*2/3)))
  layer = 1:3 # 1:3
  act.fct = c("logistic", "tanh")
  threshold = c(0.01, 0.02,0.005) #c(0.01)
  rep = c(1,5) #c(1,10)

  para = hyperaparameter(df = df, outvar = outvar, modeltype = modeltype, technique = "ann", hidden = hidden, layer = layer, act.fct =act.fct, threshold = threshold, rep =rep)
  # print(para)
  para = para[para$rep == min(para$rep),]
  para = para[1,]
  hidden = para$hidden; layer = para$layer; act.fct = para$act.fct; threshold = para$threshold; rep = para$rep
  # print(act.fct)
  fit = anntech(df = df, outvar = outvar, modeltype = modeltype, hidden = hidden, layer = layer, act.fct = act.fct, threshold = threshold, rep =rep)
  # str(fit)

  return(fit)
}
