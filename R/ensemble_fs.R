#' This function runs ensemble homogenous model
#'
#'@export
ensem_homo = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", sampletype = "bootsample", fstype = "randomfeat", maxfeatsample = ncol(traindf)-2, boots = 100, parallel = T, perfparallel = T, Trace=T, savememory = T, smart=T, ...){

  alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
              step = NA, ann = anntech, dt = NA, reg = NA, rfauto = autorftech, flasso = flasso, falasso = falasso,
              fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

  # Dataset Generation
  samplegen = c(bootsample = bootsample, randomsample = randomsample, cvsample = cvsample, allsample = allsample)
  featgen = c(allfeat = allfeat, randomfeat = randomfeat)

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
  if(Trace){print(modeltype)}

  featlist = featgen[[fstype]](df = df, outvar = outvar, maxfeatsample = maxfeatsample, boots = boots)

  traindflist = lapply(1:length(samplelist), function(x) y = list(samples = samplelist[[x]],
                                                                  features = featlist[[x]]))



  if(Trace){print("Initiate Models")}

  # Run the Model
  if(parallel){
    if(savememory){
      largeList::saveList(object = traindflist, file = "traindflist.llo", append = FALSE, compress = TRUE)
      rm(traindflist)
    }
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}

      fitlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        smartlist = future.apply::future_lapply(lb:ub, function(x){#
          # if(x>2){return(x)}
          # if(Trace){cat(x," ")}
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          newtraindf = df[trainsample, c(features,outvar)]
          newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
          fit = emb_mod(df = newtraindf, outvar = outvar, modeltype =modeltype,
                        technique = technique,...) #EnsembleFS::
        }, future.seed = T)
        # str(smartlist)
        smartlist
      })
      fitlist = purrr::flatten(fitlist)
    }
    else{
      fitlist = future.apply::future_lapply(1:length(samplelist), function(x){
                                  if(Trace){cat(x," ")}
                                  if(savememory){
                                    traindflist <- largeList::readList(file = "traindflist.llo", index = x)
                                    trainsample = traindflist[[1]]$samples
                                    features = traindflist[[1]]$features
                                  }
                                  else{
                                    trainsample = traindflist[[x]]$samples
                                    features = traindflist[[x]]$features
                                  }
                                  newtraindf = df[trainsample, c(features,outvar)]
                                  newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
                                  fit = emb_mod(df = newtraindf, outvar = outvar, modeltype =modeltype,
                                                technique = technique,...)
  }, future.seed = T)
    }
  }
  else{
    fitlist = lapply(1:length(samplelist), function(x){
      if(Trace){if(floor(x)==x){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      newtraindf = df[trainsample, c(features,outvar)]
      newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
      fit = emb_mod(df = newtraindf, outvar = outvar, modeltype =modeltype,
                    technique = technique,...)
      fit
    })#
  }

  # str(fitlist)
  # return(fitlist)
  # # Get weight
  ranklist = ranklistcreator(model_list = fitlist, techlist = technique)
  ranklist$perf = NA

  # Predict performance
  if(parallel & perfparallel){
    if(Trace){print("Initiate Performance")}
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}
      testlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        testlist = future.apply::future_sapply(lb:ub, function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
        perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                         technique = technique, model = fitlist[[x]])
        perf

      })
      })
      testlist = purrr::flatten(testlist)
      featdf = ranklist
      featdf$perf = as.numeric(testlist)

      trainlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        trainlist = future.apply::future_sapply(lb:ub, function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtraindf = df[trainsample, c(features,outvar)]
        perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                         technique = technique, model = fitlist[[x]])
        perf
      })
      })
      trainlist = purrr::flatten(trainlist)
      trainfeatdf = ranklist
      trainfeatdf$perf = as.numeric(testlist)
    }
    else{
      testlist = future.apply::future_sapply(1:length(samplelist), function(x){
      if(savememory){
        traindflist <- largeList::readList(file = "traindflist.llo", index = x)
        trainsample = traindflist[[1]]$samples
        features = traindflist[[1]]$features
      }
      else{
        trainsample = traindflist[[x]]$samples
        features = traindflist[[x]]$features
      }
      newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
      perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                       technique = technique, model = fitlist[[x]])
      perf

    })#
      featdf = ranklist
      featdf$perf = testlist

      trainlist = future.apply::future_sapply(1:length(samplelist), function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtraindf = df[trainsample, c(features,outvar)]
        perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                         technique = technique, model = fitlist[[x]])
        perf
      })#
      trainfeatdf = ranklist
      trainfeatdf$perf = trainlist
    }

  }
  else{
    testlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
                              trainsample = traindflist[[x]]$samples
                              features = traindflist[[x]]$features
                              newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
                              perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                                        technique = technique, model = fitlist[[x]])
                              perf

    })#future.apply::future_
    featdf = ranklist
    featdf$perf = testlist

    trainlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      newtraindf = df[trainsample, c(features,outvar)]
      perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                       technique = technique, model = fitlist[[x]])
      perf
    })#future.apply::future_
    trainfeatdf = ranklist
    trainfeatdf$perf = trainlist
  }

  outlist = list(fit = fitlist, df = featdf, trainperfdf = trainfeatdf, samplelist = samplelist, featlist = featlist)
  # return(fitlist)
  return(outlist)
}

#'@export
ensem_homo_int = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", sampletype = "bootsample", fstype = "randomfeat", maxfeatsample = ncol(traindf)-2, boots = 100, parallel = T, perfparallel = T, Trace=T, savememory = T, smart=T, ...){

  alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
              step = NA, ann = anntech, dt = NA, reg = NA, rfauto = autorftech, flasso = flasso, falasso = falasso,
              fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

  # Dataset Generation
  samplegen = c(bootsample = bootsample, randomsample = randomsample, cvsample = cvsample, allsample = allsample)
  featgen = c(allfeat = allfeat, randomfeat = randomfeat)

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
  if(Trace){print(modeltype)}

  featlist = featgen[[fstype]](df = df, outvar = outvar, maxfeatsample = maxfeatsample, boots = boots)

  traindflist = lapply(1:length(samplelist), function(x) y = list(samples = samplelist[[x]],
                                                                  features = featlist[[x]]))



  if(Trace){print("Initiate Models")}

  # Run the Model
  if(parallel){
    if(savememory){
      largeList::saveList(object = traindflist, file = "traindflist.llo", append = FALSE, compress = TRUE)
      rm(traindflist)
    }
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}

      fitlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        smartlist = future.apply::future_lapply(lb:ub, function(x){
          # if(x>2){return(x)}
          # if(Trace){cat(x," ")}
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          newtraindf = df[trainsample, c(features,outvar)]
          newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
          fit = AIFS2::emb_mod_int(df = newtraindf, outvar = outvar, modeltype =modeltype,
                        technique = technique,...)
        }, future.seed = T)
        # str(smartlist)
        smartlist
      })
      fitlist = purrr::flatten(fitlist)
    }
    else{
      fitlist = future.apply::future_lapply(1:length(samplelist), function(x){
        if(Trace){cat(x," ")}
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtraindf = df[trainsample, c(features,outvar)]
        newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
        fit = emb_mod_int(df = newtraindf, outvar = outvar, modeltype =modeltype,
                      technique = technique,...)
      }, future.seed = T)
    }
  }
  else{
    fitlist = lapply(1:length(samplelist), function(x){
      if(Trace){if(floor(x)==x){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      newtraindf = df[trainsample, c(features,outvar)]
      newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
      fit = emb_mod_int(df = newtraindf, outvar = outvar, modeltype =modeltype,
                    technique = technique,...)
      fit
    })#
  }

  # str(fitlist)
  # return(fitlist)
  # # Get weight
  ranklist = ranklistcreator(model_list = fitlist, techlist = technique)
  ranklist$perf = NA


  # Predict performance
  if(parallel & perfparallel){
    if(Trace){print("Initiate Performance")}
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}
      testlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        testlist = future.apply::future_sapply(lb:ub, function(x){
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
          perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                           technique = technique, model = fitlist[[x]])
          perf

        })
      })
      testlist = purrr::flatten(testlist)
      featdf = ranklist
      featdf$perf = as.numeric(testlist)

      trainlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        trainlist = future.apply::future_sapply(lb:ub, function(x){
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          newtraindf = df[trainsample, c(features,outvar)]
          perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                           technique = technique, model = fitlist[[x]])
          perf
        })
      })
      trainlist = purrr::flatten(trainlist)
      trainfeatdf = ranklist
      trainfeatdf$perf = as.numeric(testlist)
    }
    else{
      testlist = future.apply::future_sapply(1:length(samplelist), function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
        perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                         technique = technique, model = fitlist[[x]])
        perf

      })#
      featdf = ranklist
      featdf$perf = testlist

      trainlist = future.apply::future_sapply(1:length(samplelist), function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        newtraindf = df[trainsample, c(features,outvar)]
        perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                         technique = technique, model = fitlist[[x]])
        perf
      })#
      trainfeatdf = ranklist
      trainfeatdf$perf = trainlist
    }

  }
  else{
    testlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      newtestdf = df[!rownames(df) %in% trainsample, c(features,outvar)]
      perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                       technique = technique, model = fitlist[[x]])
      perf

    })#future.apply::future_
    featdf = ranklist
    featdf$perf = testlist

    trainlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      newtraindf = df[trainsample, c(features,outvar)]
      perf = perffunct_int(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                       technique = technique, model = fitlist[[x]])
      perf
    })#future.apply::future_
    trainfeatdf = ranklist
    trainfeatdf$perf = trainlist
  }

  outlist = list(fit = fitlist, df = featdf, trainperfdf = trainfeatdf, samplelist = samplelist, featlist = featlist)
  # return(fitlist)
  return(outlist)
}

#'@export
ensem_homo_cov = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", sampletype = "bootsample", fstype = "randomfeat", maxfeatsample = ncol(traindf)-2, boots = 100, parallel = T, perfparallel = T, Trace=T, savememory = T, smart=T, covs = NA, ...){

  alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
              step = NA, ann = anntech, dt = NA, reg = NA, rfauto = autorftech, flasso = flasso, falasso = falasso,
              fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

  oridf = df

  if(!is.na(covs)){df = df[,!names(df) %in% covs]}

  # Dataset Generation
  samplegen = c(bootsample = bootsample, randomsample = randomsample, cvsample = cvsample, allsample = allsample)
  featgen = c(allfeat = allfeat, randomfeat = randomfeat)

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
  if(Trace){print(modeltype)}

  featlist = featgen[[fstype]](df = df, outvar = outvar, maxfeatsample = maxfeatsample, boots = boots)

  traindflist = lapply(1:length(samplelist), function(x) y = list(samples = samplelist[[x]],
                                                                  features = featlist[[x]]))



  if(Trace){print("Initiate Models")}
  df = oridf

  # Run the Model
  if(parallel){
    if(savememory){
      largeList::saveList(object = traindflist, file = "traindflist.llo", append = FALSE, compress = TRUE)
      rm(traindflist)
    }
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}

      fitlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        smartlist = future.apply::future_lapply(lb:ub, function(x){
          # if(x>2){return(x)}
          # if(Trace){cat(x," ")}
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }

          if(!is.na(covs)){flist = c(covs,features,outvar)}
          else{flist = c(features,outvar)}

          newtraindf = df[trainsample, flist]
          newtestdf = df[!rownames(df) %in% trainsample, flist]
          fit = emb_mod_cov(df = newtraindf, outvar = outvar, modeltype =modeltype,
                        technique = technique,covs = covs,...)
        }, future.seed = T)
        # str(smartlist)
        smartlist
      })
      fitlist = purrr::flatten(fitlist)
    }
    else{
      fitlist = future.apply::future_lapply(1:length(samplelist), function(x){
        if(Trace){cat(x," ")}
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        if(!is.na(covs)){flist = c(covs,features,outvar)}
        else{flist = c(features,outvar)}

        newtraindf = df[trainsample, flist]
        newtestdf = df[!rownames(df) %in% trainsample, flist]
        fit = emb_mod_cov(df = newtraindf, outvar = outvar, modeltype =modeltype,
                      technique = technique,covs = covs,...)
      }, future.seed = T)
    }
  }
  else{
    fitlist = lapply(1:length(samplelist), function(x){
      if(Trace){if(floor(x)==x){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      if(!is.na(covs)){flist = c(covs,features,outvar)}
      else{flist = c(features,outvar)}
      # print(flist)
      # print(setdiff(flist, names(df)))
      newtraindf = df[trainsample, flist]
      newtestdf = df[!rownames(df) %in% trainsample, flist]
      # return(newtraindf)
      # print(sum(is.na(newtraindf)))
      fit = emb_mod_cov(df = newtraindf, outvar = outvar, modeltype =modeltype,
                    technique = technique, covs = covs,...)
      fit
    })#
  }

  # str(fitlist)
  # return(fitlist)
  # # Get weight
  ranklist = ranklistcreator(model_list = fitlist, techlist = technique)
  ranklist$perf = NA


  # Predict performance
  if(parallel & perfparallel){
    if(Trace){print("Initiate Performance")}
    if(smart){
      smartnumb = 10
      blocks = floor(length(samplelist)/smartnumb)
      if(length(samplelist) != blocks*smartnumb){blocks = blocks+1}
      testlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        testlist = future.apply::future_sapply(lb:ub, function(x){
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          if(!is.na(covs)){flist = c(covs,features,outvar)}
          else{flist = c(features,outvar)}

          newtestdf = df[!rownames(df) %in% trainsample, flist]
          perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                           technique = technique, model = fitlist[[x]])
          perf

        })
      })
      testlist = purrr::flatten(testlist)
      featdf = ranklist
      featdf$perf = as.numeric(testlist)

      trainlist = lapply(1:blocks, function(y){
        if(Trace){cat("Blocks",y," ")}
        lb = 1+(smartnumb*(y-1))
        ub = min(length(samplelist), smartnumb+(smartnumb*(y-1)))
        trainlist = future.apply::future_sapply(lb:ub, function(x){
          if(savememory){
            traindflist <- largeList::readList(file = "traindflist.llo", index = x)
            trainsample = traindflist[[1]]$samples
            features = traindflist[[1]]$features
          }
          else{
            trainsample = traindflist[[x]]$samples
            features = traindflist[[x]]$features
          }
          if(!is.na(covs)){flist = c(covs,features,outvar)}
          else{flist = c(features,outvar)}

          newtraindf = df[trainsample, flist]

          perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                           technique = technique, model = fitlist[[x]])
          perf
        })
      })
      trainlist = purrr::flatten(trainlist)
      trainfeatdf = ranklist
      trainfeatdf$perf = as.numeric(testlist)
    }
    else{
      testlist = future.apply::future_sapply(1:length(samplelist), function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        if(!is.na(covs)){flist = c(covs,features,outvar)}
        else{flist = c(features,outvar)}

        newtestdf = df[!rownames(df) %in% trainsample, flist]
        perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                         technique = technique, model = fitlist[[x]])
        perf

      })#
      featdf = ranklist
      featdf$perf = testlist

      trainlist = future.apply::future_sapply(1:length(samplelist), function(x){
        if(savememory){
          traindflist <- largeList::readList(file = "traindflist.llo", index = x)
          trainsample = traindflist[[1]]$samples
          features = traindflist[[1]]$features
        }
        else{
          trainsample = traindflist[[x]]$samples
          features = traindflist[[x]]$features
        }
        if(!is.na(covs)){flist = c(covs,features,outvar)}
        else{flist = c(features,outvar)}

        newtraindf = df[trainsample, flist]

        perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                         technique = technique, model = fitlist[[x]])
        perf
      })#
      trainfeatdf = ranklist
      trainfeatdf$perf = trainlist
    }

  }
  else{
    testlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      if(!is.na(covs)){flist = c(covs,features,outvar)}
      else{flist = c(features,outvar)}

      newtestdf = df[!rownames(df) %in% trainsample, flist]
      perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtestdf,
                       technique = technique, model = fitlist[[x]])
      perf

    })#future.apply::future_
    featdf = ranklist
    featdf$perf = testlist

    trainlist = sapply(1:length(samplelist), function(x){
      if(Trace){if(x/100 == floor(x/100)){cat(x," ")}}
      trainsample = traindflist[[x]]$samples
      features = traindflist[[x]]$features
      if(!is.na(covs)){flist = c(covs,features,outvar)}
      else{flist = c(features,outvar)}

      newtraindf = df[trainsample, flist]

      perf = perffunct(modeltype = modeltype, outvar = outvar, testdf = newtraindf,
                       technique = technique, model = fitlist[[x]])
      perf
    })#future.apply::future_
    trainfeatdf = ranklist
    trainfeatdf$perf = trainlist
  }

  outlist = list(fit = fitlist, df = featdf, trainperfdf = trainfeatdf, samplelist = samplelist, featlist = featlist)
  # return(fitlist)
  return(outlist)
}

