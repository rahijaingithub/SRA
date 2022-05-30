#' This function runs filter model
#'
#'@export
filter_mod = function(df = traindf, outvar = "y", modeltype = "regression", technique = "uva", ...){
  alltech = c(uva = uva)

  # Prepare Model
  fitlist = alltech[[technique]](df, outvar, modeltype)
  # fitlist = purrr::flatten(fitlist)

  # Get feature importance

  return(fitlist)
}

#'@export
uva = function(df, outvar, modeltype){

  featname = setdiff(names(df), outvar)

  if(modeltype == "regression"){
    res = future.apply::future_lapply(featname, function(x){ #

      if(which(featname == x)/10000 == floor(which(featname == x)/10000)){cat(x, " ")}

      f = as.formula(paste(outvar,"~",x))
      fit = lm(f, data = df)
    })
  }
  else if(modeltype == "logistic"){
    res = future.apply::future_lapply(featname, function(x){
      f = as.formula(paste(outvar,"~",x))
      fit = glm(f, data = df, family = "binomial")
    })
  }
  else{
    res = future.apply::future_lapply(featname, function(x){
      f = as.formula(paste("survival::Surv(SurvTime, SurvEvent)","~",x))
      fit = survival::coxph(f, data = df)
    })
  }

  return(res)
}
