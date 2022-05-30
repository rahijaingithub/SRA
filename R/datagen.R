#' Here functions generate artificial data for regression, logistic and survival
#'

inputdataset = function(varnum, setting="No_Correlation", seed=2, correlation_var=15, correlation_val=5, train_sample=500, datadist = c("random", "normal"), corr = c("fixed","random")){

  # Create the Covariance Matrix
  if(setting != "No_Correlation"){
    if(corr == "fixed"){
      Sigma=matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)
      for(i in 1:varnum){
        for(j in 1:varnum){
          if(i==j){Sigma[i,j]=10}
          else if(i<=correlation_var & j<=correlation_var & setting != "No_Correlation"){Sigma[i,j]=Sigma[j,i]=correlation_val}
          else{Sigma[i,j]=Sigma[j,i]=0}
        }
      }
      # Sigma = unlist(Sigma[lower.tri(Sigma)])
    }
    else{
      Sigma = randcorr(varnum = varnum, maxcorr = correlation_val, type = c("matrix"), seed = seed)
    }
  }
  else{
    Sigma = matrix(rep(0,varnum), nrow=varnum, ncol=varnum, byrow=F)
    diag(Sigma) <-1
    # Sigma = unlist(Sigma[lower.tri(Sigma)])
  }

  # if(max(abs(Sigma)) > 1){Sigma <- Sigma/10}

  # Generate the predictor Space
  set.seed(seed)
  ta=data.frame(MASS::mvrnorm(n = train_sample+500, rep(0, varnum), Sigma/10))

  # if(datadist == "normal"){
  #   distribution = rep("norm", varnum)
  #   data=cormatgenerator(variablesize = varnum, corr_values = Sigma, output = "both", datapoint = train_sample+500, copseed = seed, distrib = distribution)
  #   ta=data[[2]]
  # }
  # else{
  #   set.seed(seed)
  #   distribution = sample(c("beta", "norm", "unif"), varnum, replace = T)
  #   data=cormatgenerator(variablesize = varnum, corr_values = Sigma, output = "both", datapoint = train_sample+500, copseed = seed, distrib = distribution)
  #   ta=data[[2]]
  # }


  variablelist=list()
  for(i in 1:varnum){
    variablelist[[i]]=gsub(" ", "",paste("X",i))
    ta[,i]=mosaic::zscore(ta[,i])
  }
  variablelist=unlist(variablelist)
  colnames(ta) = variablelist
  # str(ta[,1:5])
  return(ta)
}

#' @export
regdataset = function(varnum, setting="No_Correlation", var=c("Both", "No_Mar", "No_Var", "No_Int"), seed=2, main_var=10, var_effect=0.5, correlation_var=15, correlation_val=5, train_sample=500, datadist = c("random", "normal"), corr = c("fixed","random")){

  # Create input parameters
  ta = inputdataset(varnum = varnum, setting=setting, seed=seed, correlation_var=correlation_var, correlation_val=correlation_val, train_sample=train_sample, datadist = datadist, corr = corr)
  ta = data.frame(ta)
  # Create the outcome Variable
  intercept=1
  betas = c(rep(var_effect,(2*main_var)-1))
  beta_value = betas[1:((2*main_var)-1)]

  if(var=="Both"){beta_value=beta_value*rep(1,(2*main_var)-1)}
  else if(var=="No_Mar"){beta_value=beta_value*c(rep(0,main_var), rep(1,main_var-1))}
  else if(var=="No_Var"){beta_value=beta_value*c(rep(0,main_var), rep(0,main_var-1))}
  else if(var=="No_int"){beta_value=beta_value*c(rep(1,main_var), rep(0,main_var-1))}
  else {beta_value = beta_value*c(rep(1,main_var), rep(0,main_var-1))}

  # Generate theinteraction term

  mar_var= paste(c(names(ta)[1:main_var]), collapse = "+")
  int_var= paste(names(ta)[1:(main_var-1)],"*" , names(ta)[2:main_var], collapse = " + ")

  f = as.formula(paste("~",mar_var ," +" , int_var))
  # print(f)
  main_mat = model.matrix(f, ta)
  # print(main_mat[1:5,])
  # Get Outcome
  set.seed(2)
  random_value = rnorm(n=train_sample+500, mean=0, sd=0.25)
  # coef_value = apply(main_mat,1, function(x) {sum(x*c(intercept, beta_value))})

  ta$y = main_mat %*% c(intercept, beta_value) + random_value
  # ta$y  = coef_value + random_value

  set.seed(2)
  index=sample(1:nrow(ta), train_sample, replace = F)
  traindf=ta[index,]
  validationdf=ta[-index,]
  return(list(train=traindf,test=validationdf, fulldf = ta))
  }

#' @export
bindataset = function(varnum, setting="No_Correlation", var=c("Both", "No_Mar", "No_Var", "No_Int"), seed=2, main_var=10, var_effect=0.5, correlation_var=15, correlation_val=5, train_sample=500, datadist = c("random", "normal"), corr = c("fixed","random"), outcome = c("balanced", 0.3)){

  # Create input parameters
  ta = inputdataset(varnum = varnum, setting=setting, seed=seed, correlation_var=correlation_var, correlation_val=correlation_val, train_sample=9500, datadist = datadist, corr = corr)
  ta = data.frame(ta)

  # Create the outcome Variable
  intercept=1
  betas = c(rep(var_effect,(2*main_var)-1))
  beta_value = betas[1:((2*main_var)-1)]

  if(var=="Both"){beta_value=beta_value*rep(1,(2*main_var)-1)}
  else if(var=="No_Mar"){beta_value=beta_value*c(rep(0,main_var), rep(1,main_var-1))}
  else if(var=="No_Var"){beta_value=beta_value*c(rep(0,main_var), rep(0,main_var-1))}
  else if(var=="No_int"){beta_value=beta_value*c(rep(1,main_var), rep(0,main_var-1))}
  else {beta_value = beta_value*c(rep(1,main_var), rep(0,main_var-1))}

  # Generate theinteraction term

  mar_var= paste(c(names(ta)[1:main_var]), collapse = "+")
  int_var= paste(names(ta)[1:(main_var-1)],"*" , names(ta)[2:main_var], collapse = " + ")

  f = as.formula(paste("~",mar_var ," +" , int_var))
  main_mat = model.matrix(f, ta)

  # Get Outcome
  set.seed(2)
  random_value = rnorm(n=9500+500, mean=0, sd=0.25)
  ta$y = main_mat %*% c(intercept, beta_value) + random_value

  pr = 1/(1+exp(-ta$y)) # get probability of success
  ta$y = rbinom(10000,1,pr)      # bernoulli response variable

  ta$y = as.factor(ta$y)
  successrow = rownames(ta[ta$y == 1,])
  failurerow = rownames(ta[ta$y == 0,])

  set.seed(2)
  if(outcome == "balanced"){
    success_row = sample(successrow,(train_sample+500)*0.5, replace = F)
    failure_row = sample(failurerow,(train_sample+500)*0.5, replace = F)
    traindf = ta[c(success_row[1:(train_sample*0.5)], failure_row[1:(train_sample*0.5)]),]
    validationdf = ta[c(success_row[-1:-(train_sample*0.5)], failure_row[-1:-(train_sample*0.5)]),]
  }
  else{
    success_row = sample(successrow,(train_sample+500)*outcome, replace = F)
    failure_row = sample(failurerow,(train_sample+500)*outcome, replace = F)
    traindf = ta[c(success_row[1:(train_sample*outcome)], failure_row[1:(train_sample*outcome)]),]
    validationdf = ta[c(success_row[-1:-(train_sample*outcome)], failure_row[-1:-(train_sample*outcome)]),]
  }

  # str(ta)
  return(list(train=traindf,test=validationdf, fulldf = ta))
}

#' @export
survdataset = function(varnum, setting="No_Correlation", var=c("Both", "No_Mar", "No_Var", "No_Int"), seed=2, main_var=10, var_effect=0.5, correlation_var=15, correlation_val=5, train_sample=500, datadist = c("random", "normal"), corr = c("fixed","random"), lambdas= 0.6 , gammas = 1.5, maxt=NULL, censoringrate=0.7, censortime = 0.5, model_distribution="exponential"){

  # Create input parameters
  ta = inputdataset(varnum = varnum, setting=setting, seed=seed, correlation_var=correlation_var, correlation_val=correlation_val, train_sample=train_sample, datadist = datadist, corr = corr)
  ta = data.frame(ta)

  # Create the outcome Variable
  intercept=1
  betas = c(rep(var_effect,(2*main_var)-1))
  beta_value = betas[1:((2*main_var)-1)]

  if(var=="Both"){beta_value=beta_value*rep(1,(2*main_var)-1)}
  else if(var=="No_Mar"){beta_value=beta_value*c(rep(0,main_var), rep(1,main_var-1))}
  else if(var=="No_Var"){beta_value=beta_value*c(rep(0,main_var), rep(0,main_var-1))}
  else if(var=="No_int"){beta_value=beta_value*c(rep(1,main_var), rep(0,main_var-1))}
  else {beta_value = beta_value*c(rep(1,main_var), rep(0,main_var-1))}

  # Generate theinteraction term

  mar_var= paste(c(names(ta)[1:main_var]), collapse = "+")
  int_var= paste(names(ta)[1:(main_var-1)],"*" , names(ta)[2:main_var], collapse = " + ")

  f = as.formula(paste("~",mar_var ," +" , int_var))
  # print(f)
  main_mat = model.matrix(f, ta)
  # print(main_mat[1:5,])

  # Prepare the outcome
  covs = data.frame(main_mat[,-1], stringsAsFactors = F)
  inputbeta = beta_value
  names(inputbeta) = names(covs)

  set.seed(2)
  #covs <- data.frame(id = 1:100, trt = stats::rnorm(100, 1, 0.5))
  if(model_distribution=="exponential"){
    s1 <- simsurv::simsurv(dist = model_distribution, lambdas = lambdas, betas = inputbeta, x = covs, maxt = maxt)}
  else{s1 <- simsurv::simsurv(dist = model_distribution, lambdas = lambdas, gammas=gammas, betas = inputbeta, x = covs, maxt = maxt)}

  # Randomly censor the rows
  set.seed(2)
  censor_rows = sample(1:nrow(s1), floor(nrow(s1)*censoringrate), replace = F)
  s1$eventtime[censor_rows] = s1$eventtime[censor_rows]*censortime
  s1$status[censor_rows] = 0
  # plot(survfit(Surv(eventtime, status)~1, data = s1))
  # print(sum(s1$status)/nrow(s1))

  dataset= cbind(ta,s1)
  #str(dataset)
  dataset$id= NULL
  names(dataset)[names(dataset) == "eventtime"] = "SurvTime"
  names(dataset)[names(dataset) == "status"] = "SurvEvent"
  #str(dataset)

  ## Create stratified Test and Training data
  train_rows = train_sample/(train_sample + 500)
  set.seed(2)
  SURV_1=dataset[which(dataset$SurvEvent==1),]
  S_1=sample(rownames(dataset[which(dataset$SurvEvent==1),]), floor(train_rows*nrow(SURV_1)), replace=F)
  SURV_0=dataset[which(dataset$SurvEvent==0),]
  S_0=sample(rownames(dataset[which(dataset$SurvEvent==0),]), floor(train_rows*nrow(SURV_0)), replace=F)
  index <- as.numeric(c(S_1,S_0))


  traindf=dataset[index,]
  validationdf=dataset[-index,]

  return(list(train=traindf,test=validationdf, fulldf = ta))
}

## Remove the variables with missing data
missing_data=function(inputdf, missingpercent){
  # Determine the missing data distribution across the dataset
  Missing_data=sapply(inputdf, function(x) 100*sum(is.na(x))/nrow(inputdf))
  # Remove variables with more than "missingpercent"
  var=attributes(Missing_data)$names[Missing_data>missingpercent]
  inputdf[,var]=NULL
  return(inputdf)
}
## Remove the variables with zero or NA sd
sdrem = function(inputdf, minsd = 0){

  df = apply(inputdf, 2, function(x){ y = sd(x, na.rm = T); if(!is.na(y) & y>minsd){x}else{NULL}})
  if(is.list(df)){inputdf = do.call(cbind, df)}
  else{inputdf = df}
  inputdf = data.frame(inputdf)

  return(inputdf)
}
## Remove the variables with low cv
cvrem = function(inputdf, mincv = 0.1){

  df = apply(inputdf, 2, function(x){ y = sd(x, na.rm = T); my = mean(x, na.rm = T); if(!is.na(y/my) & y/my>mincv){x}else{NULL}})
  if(is.list(df)){inputdf = do.call(cbind, df)}
  else{inputdf = df}
  inputdf = data.frame(inputdf)

  return(inputdf)
}
## Remove the variables with low variation
lowvar=function(inputdf, minfreq=100/10){
  al= unlist(apply(inputdf,2,function(x) caret::nearZeroVar(x, freqCut = minfreq, saveMetrics = T)[4]))
  drop_var_2=names(al)[which(al==TRUE)]
  drop_var_2 = stringr::str_replace(drop_var_2,".nzv","")
  newdf=inputdf[,setdiff(names(inputdf),drop_var_2)]
  newdf = data.frame(newdf)
  names(newdf) = setdiff(names(inputdf),drop_var_2)
  return(newdf)
}

## Scale the data
datascale = function(inputdf = df, scaletype = c("z", "ranknormal"), outvar = names(df)[ncol(df)], modeltype = "regression"){
  if(modeltype == "regression"){featdf = inputdf}
  else{featdf = inputdf[,!names(inputdf) %in% outvar]}

  if(scaletype == "z"){featdf = data.frame(scale(featdf))}
  else{
    featdf = apply(featf,2,function(x) bestNormalize::orderNorm(x)$x.t)
  }

  if(modeltype == "regression"){df = featdf}
  else{
    df = cbind(featdf, inputdf[,outvar])
    names(df)[ncol(df)] = outvar
  }

  # str(df)
  return(df)
}
## Reduce it to maximum columns
col_reduce =  function(inputdf, outvar, reducetech = c("uva", "lasso"), modeltype = "regression", rankcutoff = 100){

  if(reducetech == "lasso"){
    res = emb_mod(df = inputdf, outvar = outvar, modeltype = modeltype, technique = "lasso")
    finrank=loopdf(df = res, tech = "lasso")
    selfeat = finrank$variable[finrank$importance != 0]
  }
  else{
    res = filter_mod(df = inputdf, outvar = outvar, modeltype = modeltype, technique = reducetech)
    featdf = lapply(res, function(x){#future.apply::future_
                                     sumfit = summary(x)
                                     featcoef = sumfit$coefficients
                                     pvalcol = grepl("Pr",colnames(featcoef))
                                     rownumb = which(rownames(featcoef) != "(Intercept)")
                                     featname = rownames(featcoef)[rownumb]
                                     pval = featcoef[rownumb,pvalcol]
                                     imp = log10(1/pval)

                                     df = c(variable= featname, importance = imp)
                                    })
    finrank = data.frame(do.call(rbind, featdf), stringsAsFactors = F)
    finrank = finrank[order(finrank$importance, decreasing = T),]
    selfeat = finrank$variable[1:rankcutoff]
  }

  inputdf = inputdf[,c(selfeat, outvar)]
  return(inputdf)
}


#' @export
realdataprep = function(filename = NA, df = NA, miss_percent = 90, scaletype = c("z", "ranknormal"), minfreq=100/10, maxcol = 100, cv_y_miss = 5, mincol = 15, modeltype = "regression", mincv =0.1, reducetech = "uva"){
  print(filename)
  # Get df
  if(is.na(df)){
    df = read.csv(filename)
    names(df) = stringi::stri_replace_all_fixed(names(df), ".","")
    }

  # Remove duplicated columns
  df = df[!duplicated(as.list(df))]

  # Extra CNV processing
  if(!is.na(filename)){
    if(grepl("CNV", filename)){
      # print(filename)
      df = df[,-1]

      # Remove Columns with zero sd
      df =  sdrem(inputdf = df, minsd = 0)

      # Remove columns with missing values
      df = missing_data(inputdf = df, missingpercent = cv_y_miss)

      # Remove columns with low coefficient of variation
      df =  lowvar(inputdf = df, minfreq=minfreq)

      # List outcome columns
      outnames = names(df)[grepl("GDSC", names(df))]
      print(length(outnames))

      # df = data.matrix(df)
      df = data.table::data.table(df)

      # Create df with each outcome
      dflist = lapply(outnames[c(10, 41, 45, 76, 97, 108, 111, 114, 115, 121, 134, 157, 166, 178, 198)], function(x){                                                cat(which(outnames == x), " ")
                                             # Create df with each outcome
                                             y = setdiff(colnames(df),outnames);
                                             # z = df[,c(y,x)]
                                             out = df[,x]
                                             out = which(!is.na(out))
                                             # Create df with complete data of outcome
                                             z = df[out,c(y,x)]

                                             # Remove columns with missing values
                                             z = missing_data(inputdf = z, missingpercent = miss_percent)

                                             # Create df with complete data
                                             z = z[complete.cases(z),];

                                             if(nrow(z) < 2 | ncol(z) <2 | !all(x %in% names(z))){return(NULL)}

                                             # Scale the data
                                             z = datascale(inputdf = z, scaletype = "z",
                                                            outvar = x,
                                                            modeltype = modeltype)

                                             # Remove duplicated columns
                                             z = z[!duplicated(as.list(z))]

                                             # Select data upto max columns
                                             if(ncol(z) > maxcol){
                                               z[,ncol(z)] = bestNormalize::orderNorm(z[,ncol(z)])$x.t
                                               z= data.frame(z)
                                               z = col_reduce(inputdf = z, outvar = x, reducetech = reducetech,
                                                              modeltype = modeltype, rankcutoff = maxcol)}
                                             else{z}

                                             # Select data with more than minimum columns
                                             if(ncol(z) >= mincol){z}else{z=NULL}
                                             cat(dim(z), " ");
                                             # Find the most relevant Drug
                                             if(ncol(z) >50){
                                               z[,ncol(z)] = bestNormalize::orderNorm(z[,ncol(z)])$x.t
                                               z= data.frame(z)

                                               drdf = z[,-ncol(z)]

                                               # # Perform PCA
                                               # pca = prcomp(drdf, scale. = F, center = F, rank. = 2)
                                               # pca_corr = apply(pca$x,2,function(x) (cor(x,z[,ncol(z)]))^2)
                                               # cat(pca_corr, " ")

                                               # # Perform UMAP
                                               # uma = umap::umap(drdf, n_components = 2)
                                               # umap_corr = apply(uma$layout,2,function(x) (cor(x,z[,ncol(z)]))^2)
                                               # cat(umap_corr, " ")

                                               # Perform PLS
                                               plsreg = plsdepot::plsreg1(predictors = drdf, response = z[,ncol(z)],
                                                                          comps = 2, crosval = T)
                                               plsreg_corr = apply(plsreg$x.scores,2,function(x) (cor(x,z[,ncol(z)]))^2)
                                               # cat(plsreg_corr, " ")

                                               corlist = c(plsreg_corr) #pca_corr,umap_corr,
                                               if(max(corlist) < 0.175){Z = NULL}
                                               else{ z}#print(max(corlist));
                                               # plot(corlist)
                                             }
                                             z
                                             })

      # Top Priority
      toppri = c(3,15)

      res_df = lapply(dflist[toppri], function(mdf){
          plsreg = plsdepot::plsreg1(predictors = mdf[,-ncol(mdf)], response = mdf[,ncol(mdf)], comps = 2, crosval = T)
          plsreg_corr = apply(plsreg$x.scores,2,function(x) (cor(x,mdf[,ncol(mdf)]))^2)
          # print(plsreg_corr)

          # Get the features importance sequence
          featseq = names(sort(abs(plsreg$x.loads[,1])))
          outvar = names(mdf)[ncol(mdf)]
          newdf = mdf[,c(featseq[-1:-(148*100)], outvar)]
          # str(newdf)
          newdf

        })

      df = res_df
    }
    else{
      if(modeltype != "survival"){outvar = names(df)[ncol(df)]}
      else{outvar = names(df)[(ncol(df)-1):ncol(df)]}

      df = df[complete.cases(df[,outvar]),]
      outdf = data.frame(df[,outvar])
      names(outdf) = outvar


      # Remove Columns with zero sd
      df =  sdrem(inputdf = df, minsd = 0)

      if(!all(outvar %in% names(df))){return(NULL)}

      # Remove columns with missing values
      df = missing_data(inputdf = df, missingpercent = miss_percent)
      if(ncol(df) == length(outvar)){return(NULL)}

      # Remove columns with low variance
      inputdf = data.frame(df[,!names(df) %in% outvar])
      names(inputdf) = setdiff(names(df), outvar)
      df = inputdf
      df = lowvar(inputdf = df, minfreq=minfreq)
      if(ncol(df) == 0){return(NULL)}
      df = cbind(df, outdf)

      # Get complete cases
      df = df[complete.cases(df),]

      if(length(df) == 0){return(NULL)}
      if(is.na(df)){return(NULL)}
      if(is.null(dim(df))){return(NULL)}
      if(nrow(df) < 2 | ncol(df) <2 | !all(outvar %in% names(df))){return(NULL)}

      # Scale the data
      df = datascale(inputdf = df, scaletype = "z", outvar = names(df)[ncol(df)], modeltype = modeltype)

      # Remove duplicated columns
      df = df[!duplicated(as.list(df))]

      # Select data upto max columns
      if(ncol(df) > maxcol){df = col_reduce(inputdf = df, outvar = names(df)[ncol(df)], reducetech = reducetech, modeltype = modeltype, rankcutoff = maxcol)}
      else{df}

      # Select data with more than minimum columns
      if(ncol(df) >= mincol){df}else{df=NULL}
    }
  }
  else{
    if(modeltype != "survival"){outvar = names(df)[ncol(df)]}
    else{outvar = names(df)[(ncol(df)-1):ncol(df)]}

    df = df[complete.cases(df[,outvar]),]
    outdf = data.frame(df[,outvar])
    names(outdf) = outvar


    # Remove Columns with zero sd
    df =  sdrem(inputdf = df, minsd = 0)

    if(!all(outvar %in% names(df))){return(NULL)}

    # Remove columns with missing values
    df = missing_data(inputdf = df, missingpercent = miss_percent)
    if(ncol(df) == length(outvar)){return(NULL)}

    # Remove columns with low variance
    inputdf = data.frame(df[,!names(df) %in% outvar])
    names(inputdf) = setdiff(names(df), outvar)
    df = inputdf
    df = lowvar(inputdf = df, minfreq=minfreq)
    if(ncol(df) == 0){return(NULL)}
    df = cbind(df, outdf)

    # Get complete cases
    df = df[complete.cases(df),]

    if(length(df) == 0){return(NULL)}
    if(is.na(df)){return(NULL)}
    if(is.null(dim(df))){return(NULL)}
    if(nrow(df) < 2 | ncol(df) <2 | !all(outvar %in% names(df))){return(NULL)}

    # Scale the data
    df = datascale(inputdf = df, scaletype = "z", outvar = names(df)[ncol(df)], modeltype = modeltype)

    # Remove duplicated columns
    df = df[!duplicated(as.list(df))]

    # Select data upto max columns
    if(ncol(df) > maxcol){df = col_reduce(inputdf = df, outvar = names(df)[ncol(df)], reducetech = reducetech, modeltype = modeltype, rankcutoff = maxcol)}
    else{df}

    # Select data with more than minimum columns
    if(ncol(df) >= mincol){df}else{df=NULL}

  }


  # finaldf
  return(df)
}

#' @export
realdataset = function(traindf = NA, testdf = NA, outvar, train_testratio = 0.8, modeltype, removevar, seed = 1){
  if(is.na(testdf)){
    # str(traindf)
    # Remove unwanted features and clean the df
    df = realdatacleaner(df = traindf, removevar = removevar, outvar=outvar, modeltype = modeltype)

    # Split Training data
    set.seed(seed)
    sam = sample(rownames(df), floor(nrow(df)*train_testratio), replace = F)

    traindf = df[rownames(df) %in% sam,]
    testdf = df[!rownames(df) %in% sam,]
  }
  else{
    # Remove unwanted features and clean the df
    traindf = realdatacleaner(df = traindf, removevar = removevar, outvar=outvar, modeltype = modeltype)
    testdf = realdatacleaner(df = testdf, removevar = removevar, outvar=outvar, modeltype = modeltype)
  }

  return(list(traindf = traindf, testdf = testdf))
}



# Support Functions
randcorr = function(varnum = 5, maxcorr = 0.5, type = c("upperlist", "matrix"), seed = 1){
  a = varnum
  set.seed(seed)
  cormat <- matrix(runif(a*a, min=-maxcorr,max=maxcorr),a,a)
  diag(cormat) <- 1
  cormat[lower.tri(cormat)] = t(cormat)[lower.tri(cormat)]

  if(type == "upperlist"){
    cormat = unlist(cormat[lower.tri(cormat)])
  }

  return(cormat)
}
param_dist = function(varnum = 5, distrib = NA, ...){

  # Default Feature Values
  fullvarlist = list(shape1=7, shape2=2, min = 0, max =2, mean = 0.1, sd = 1)
  varlist= list(...)
  ## Replace standard variable values with user defined values
  if(length(varlist) >0){
    targetvar = names(fullvarlist)[names(fullvarlist) %in% names(varlist)]
    dumpvar = sapply(targetvar, function(x) fullvarlist[[x]]<<-varlist[[x]])
  }

  ## All Standard Features
  shape1= fullvarlist$shape1;
  shape2=fullvarlist$shape2;
  min = fullvarlist$min;
  max = fullvarlist$max;
  mean = fullvarlist$mean;
  sd = fullvarlist$sd

  # Prepare distribution

  if(is.na(distrib)){distrib=c("beta", "norm", "unif", rep("norm",max(varnum-3,0)))}

  distrib_para=lapply(distrib, function(x) {
    seeder = sample(seq(1,5*varnum),1)
    set.seed(seeder);
    if (x=="beta"){out = list(shape1= shape1, shape2=shape2)}
    else if(x=="unif"){out = list(min= min, max= max)}
    else {
      out = list(mean=mean, sd=sd)}
    out
  })

  return(list(distrib, distrib_para))
}

cormatgenerator=function(variablesize=non_zerocorr, corr_values=corr_list, output=c("cormat", "data", "both"),
                         datapoint=1000, copseed=1, ...){

  # Default Feature Values
  fullvarlist = list(shape1=7, shape2=2, min = 0, max =2, mean = 0.1, sd = 1, distrib = NA)
  varlist= list(...)
  ## Replace standard variable values with user defined values
  if(length(varlist) >0){
    targetvar = names(fullvarlist)[names(fullvarlist) %in% names(varlist)]
    dumpvar = sapply(targetvar, function(x) fullvarlist[[x]]<<-varlist[[x]])
  }

  ## All Standard Features
  shape1= fullvarlist$shape1;
  shape2=fullvarlist$shape2;
  min = fullvarlist$min;
  max = fullvarlist$max;
  mean = fullvarlist$mean;
  sd = fullvarlist$sd
  distrib = fullvarlist$distrib

  # print(c(variablesize,corr_values))
  # if(!is.null(copseed)==T){set.seed(copseed)},
  set.seed(copseed)
  myCop=copula::normalCopula(param=corr_values,dim = variablesize, dispstr = "un")
  distribution_para = param_dist(varnum = variablesize, distrib = distrib, shape1 = shape1, shape2 = shape2, min = min, max = max, mean = mean, sd = sd)
  distrib = distribution_para[[1]]
  distrib_para = distribution_para[[2]]
  # print(distrib_para)
  myMvd <- copula::mvdc(copula=myCop, margins=distrib, paramMargins=distrib_para )
  #set.seed(1)
  Z2 <- copula::rMvdc(datapoint,myMvd)
  cor_mat=cor(Z2)
  if(output=="cormat"){return(cor_mat)}
  else if(output=="data"){return(Z2)}
  else{return(list(cor_mat,Z2))}
}

realdatacleaner = function(df, removevar = NA, outvar="y", modeltype = "logistic"){

  # remove unwanted columns and reorganise the outcome feature
  df = df[, !names(df) %in% removevar]
  if(length(outvar) > 1){
    dffeat = df[,!names(df) %in% outvar]
    dfout = df[,outvar]
    names(dfout) = c("SurvTime", "SurvEvent")
    df = cbind(dffeat, dfout)
  }
  else{
    dffeat = df[,!names(df) %in% outvar]
    dfout = data.frame(df[,outvar], stringsAsFactors = F)
    names(dfout) = "y"
    if(modeltype == "logistic"){dfout$y = as.factor(dfout$y)}
    df = cbind(dffeat, dfout)
  }
  return(df)
}


# Missing data generator
missgen = function(df, outvar, missdata = 0.8, type = c("scr", "ncr"), scr = 50){

  if("scr" %in% type){
    completedf = df[1:scr,]
    missdf = df[-1:-scr,]
  }
  else{missdf = df}

  # Columns containing outcome variable names
  outvarcol = which(names(df) %in% outvar)

  # Get inputdf
  inputdf = missdf[,!names(missdf) %in% outvar]

  # Create missing values in all rows
  inputdflist = lapply(1:nrow(inputdf), function(x){
                                          set.seed(x);
                                          cindex = sample(1:ncol(inputdf), floor(missdata*ncol(inputdf)), replace = F)
                                          rowval = inputdf[x,]
                                          rowval[cindex] = NA
                                          rowval
                                          })
  inputdf = do.call(rbind, inputdflist)

  # Combine outvar with inputdf
  missdf = cbind(inputdf,missdf[,names(df) %in% outvar])

  names(missdf)[outvarcol] = outvar
  # Create final df with complete and missdf
  if("scr" %in% type){finaldf = rbind(completedf, missdf)}
  else{finaldf = missdf}

  return(finaldf)
}
