#' Hyperparameter optimization functions
#'
#' These function optimize the different techniques for Survival, Binary or Continuous
#'
#' @export
hyperaparameter = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", ...){
  alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
              step = NA, ann = anntech, dt = NA, reg= NA)
  # Define the parameters
  ## Split the data
  samplegen = c(bootsample = bootsample, randomsample = randomsample, cvsample = cvsample, allsample = allsample)
  sampletype = "randomsample"
  boots = 1
  if(modeltype == "regression"){
    samplelist = samplegen[[sampletype]](df = df, boots = boots)
  }
  else{
    stratvar = ifelse(length(outvar) ==1, outvar, "SurvEvent")
    df0 = df[df[,stratvar] == 0,]
    samplelist0 = samplegen[[sampletype]](df = df0, boots = boots)
    df1 = df[df[,stratvar] == 1,]
    samplelist1 = samplegen[[sampletype]](df = df1, boots = boots)
    samplelist = lapply(1:length(samplelist0), function(x) c(samplelist0[[x]], samplelist1[[x]]))
  }
  # str(samplelist)
  trainsample = samplelist[[1]]
  traindf = df[trainsample,]
  testdf = df[!rownames(df) %in% trainsample,]

  ## Check the CV
  arglist = list(...)
  # print(arglist)
  ## Check the Technique Parameters
  hyperdf = expand.grid(arglist, stringsAsFactors = F)
  # str(hyperdf)

  if(any(technique %in% c("rf"))){
    perf = sapply(1:nrow(hyperdf), function(x){
            # if(x/100 == floor(x/100)){cat(x, " "); print(hyperdf[x,])}
            fit = tryCatch({do.call(alltech[[technique]],c(hyperdf[x,], list(df = traindf, outvar=outvar, modeltype = modeltype)))}, error = function(cond) return(NA))
            if(is.na(fit)){perfmet = NA}
            else{perfmet = perffunct(modeltype =  modeltype, outvar = outvar, testdf = testdf, pred = fit$predicted, technique = technique, model = NA)}
            perfmet
    })
    parameters = hyperdf[which(perf == max(perf, na.rm = T)),]
    parameters = parameters[1,]
  }
  else if(any(technique %in% c("ann"))){
    perf = sapply(1:nrow(hyperdf), function(x){
      # print(hyperdf[x,])
      fit = tryCatch({do.call(alltech[[technique]],c(hyperdf[x,], list(df = traindf, outvar=outvar, modeltype = modeltype)))}, error = function(cond) return(NA))

      if(is.na(fit)){perfmet = NA}
      else{
        if(is.null(fit$net.result)){perfmet = NA}
        else{
          Predict_try = neuralnet::compute(fit,testdf)
          pred_try <- Predict_try$net.result

          perfmet = perffunct(modeltype =  modeltype, outvar = outvar, testdf = testdf, pred = pred_try, technique = technique, model = NA)
          # perfmet = sapply(1:length(fit$net.result), function (i){
          #   if(is.null(fit$net.result[[i]])){perfmet = NA}
          #   else{}
          #   perfmet
          # })
          # perfmet = max(perfmet, na.rm = T)
          perfmet
        }
      }
      perfmet
    })
    parameters = hyperdf[which(perf == max(perf, na.rm = T)),]
  }
  else if (technique %in% c("step", "dt")){fit = NA}
  else{
    if(!names(arglist) %in% "cv"){hyperdf$cv = 5}
    fit = alltech[[technique]](df, outvar, modeltype, ...)
  }
  # print(parameters)
  # fit = do.call(alltech[[technique]],c(parameters, list(df = df, outvar=outvar, modeltype = modeltype)))
  # print(fit)
  return(parameters)
}

