#' This function runs rank aggregation
#'
#'@export
rankaggfun = function(rankmethod, model_list = NA, techlist = NA, ranklist = NA, covs = NA, ...){
  # Rank Aggregation
  rankagg = c(rulebase = rulebase, unsuprankagg = unsuprankagg, suprankagg = suprankagg)

  if(!all(is.na(model_list))){
    ranklist = ranklistcreator(model_list = model_list, techlist = techlist)
  }

  if(!is.na(covs)){ranklist = ranklist[,!names(ranklist) %in% covs]}

  # Feature Ranking
  finrank = rankagg[[rankmethod]](ranklist, ...)

  return(finrank)
}

ranklistcreator = function(model_list = fitlist, techlist = NA){
  techlist = unlist(techlist)

  if(is.na(techlist)){print("Cannot Process, technique argument is not defined"); res = NA}
  else if(length(techlist) >1 & length(techlist) < length(model_list)){print("Cannot Process, technique argument is length is more than 1 and less than length of model_list"); res = NA}
  else{
    if(length(techlist) == 1){techlist = rep(techlist, length(model_list))}
    ranklist = lapply(1:length(model_list), function(x){
      # cat(x, " ")
      # str(model_list[[x]])
                                            featrank = loopdf(df = model_list[[x]], tech = techlist[x])
                                            names(featrank) = c("variable", paste("normwt",x,sep = ""))
                                            # featrank$rank = rank(featrank$normwt)
                                            # names(featrank)[3] = c(paste("rank",x,sep = ""))
                                            featrank
                                            })
    res = ranklist
    res = Reduce(function(...) merge(..., by="variable", all=TRUE), res)
    tres = t(res[,-1])
    colnames(tres) = res[,1]
    tres = data.frame(tres)
    res = tres
  }

  resrank2 = res
  names(resrank2) = var_organise(names(res), ".")
  varlist = table(var_organise(names(res), "."))

  trial = lapply(1:length(varlist), function(x) {
    if(varlist[x] == 1){
      varname = names(varlist)[x]
      y = data.frame(resrank2[,varname])
      names(y) = varname
      y
    }
    else{
      varname = names(varlist)[x]
      meany = rowMeans(resrank2[,grepl(varname, names(resrank2))], na.rm = T)
      meany[is.nan(meany)] = NA
      y = data.frame(meany)
      names(y) = varname
      y
    }
  })
  res = do.call(cbind, trial)
  # str(res)
  return(res)
}

#'@export
rulebase = function(ranklist, rulemethod = c("min", "max", "mean", "median", "freq", "sd", "coefv", "ttest", "wmw"), perfvar = NA){
  if(!is.na(perfvar)){ranklist = ranklist[ ,!names(ranklist) %in% perfvar]}

  ruletech = c(min = min, max = max, mean = mean, median = median, freq = freq, sd = sdev, coefv = coefval, ttest = ttest, wmw = wtest, rra = rra, stuart = stuart, rankmean = rankmean, rankmedian = rankmedian, rankgeomean = rankgeomean, rankmin = rankmin)

  if(rulemethod %in% c("min", "max", "mean", "median", "freq", "sd", "coefv", "ttest", "wmw")){
    featwt = sapply(1:ncol(ranklist), function(x){y = ranklist[,x]; ruletech[[rulemethod]](y, na.rm = T)})
    featwt = abs(featwt)
    featrank = data.frame(variable = names(ranklist), imp = as.numeric(featwt), stringsAsFactors = F)
    featrank$rank = rank(1/featrank$imp)
  }
  else{
    featrank = ruletech[[rulemethod]](y=ranklist, na.rm = T)
    # print("Before Log")
    # print(featrank$imp)
    # featrank$imp = log10(featrank$imp)
    featrank$imp = asinh(featrank$imp)
  }

  return(featrank)
}

freq = function(y, na.rm = T){
  if(na.rm){
    y = y[!is.na(y)]
  }
  n = length(y)

  varfreq = length(y[y!=0])/n
  # cat(n, " ", varfreq," ")
  return(varfreq)
}
ttest = function(y, na.rm  = T){
  # print(unique(y)[!is.na(unique(y))])

  z = unique(y)
  if(length(z[!is.na(z)]) <= 1){y = sapply(y,function(x) x + sample(c(1e-9, -1e-9),1))}
  if(length(unique(y[!is.na(y)])) <= 1){
    y = sapply(y,function(x) x + sample(c(-9:-1, 1:9)*1e-9,1, replace = T))
  }
  if(length(unique(y[!is.na(y)])) <= 1){return(0)}
  tval = t.test(y)$statistic

  return(as.numeric(tval))
}
wtest = function(y, na.rm  = T){
  wval = wilcox.test(y)$statistic
  return(as.numeric(wval))
}
coefval = function(y, na.rm  = T){
  m = mean(y, na.rm  = T)
  s = sd(y, na.rm  = T)
  cv = s/m
  cv = 1/cv
  return(as.numeric(cv))
}
sdev = function(y, na.rm  = T){
  # print(T)
  s = sd(y, na.rm  = T)
  s = 1/s
  return(as.numeric(s))
}

rrapack = function(y, rankprep = T, ...){

  # Check Input structure
  checkval = as.numeric(unlist(y))
  if(min(checkval, na.rm = T) <0 | max(checkval, na.rm = T)>1){normy = y}else{normy = NA}

  if(any(class(y) %in% "matrix") & is.na(normy)){r = y}else{r = NA}

  # Prepare Rank and/or matrix
  if(is.na(r) | rankprep == T){
    # Prepare the rank: Each row is a sample, Largest value gets highest rank
    rows = nrow(y)
    ranklist = lapply(1:rows, function(x){y[x,] = 1/abs(y[x,]); rankorder = rank(y[x,], na.last = NA); featrank = names(sort(rankorder))})
    r = RobustRankAggreg::rankMatrix(ranklist)
  }

  rra = RobustRankAggreg::aggregateRanks(rmat = r, ...)

  return(rra)
}
rra = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "RRA")
  featrank = data.frame(variable = res[,1], imp = 1/as.numeric(res[,2]), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
stuart = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "stuart")
  featrank = data.frame(variable = res[,1], imp = 1/abs(as.numeric(res[,2])), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
rankmean = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "mean")
  featrank = data.frame(variable = res[,1], imp = 1/as.numeric(res[,2]), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
rankmedian = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "median")
  featrank = data.frame(variable = res[,1], imp = 1/as.numeric(res[,2]), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
rankgeomean = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "geom.mean")
  featrank = data.frame(variable = res[,1], imp = 1/as.numeric(res[,2]), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
rankmin = function(y, na.rm = T){
  res = rrapack(y = y, rankprep = T, method = "min")
  featrank = data.frame(variable = res[,1], imp = 1/as.numeric(res[,2]), stringsAsFactors = F)
  featrank$rank = rank(1/featrank$imp)
  return(featrank)
}
km = function(y, na.rm = T){

  y[is.na(y)] = 0
  # Remove Performance column ("Perf") if any
  y[,"perf"] = NULL

  featname = names(y)
  # Absolute value
  y = abs(y)

  # Transpose
  y = data.frame(t(y))

  # Max Scaling
  y = purrr::reduce(purrr::map(y, function(x) {x/max(x)}), data.frame)
  y[is.na(y)] = 0;

  if(class(y) != "data.frame"){y = data.frame(val = y)}
  rownames(y) = featname

  # str(y)
  if(ncol(y) >1){
    if(nrow(y[!duplicated(y),]) > 2 | (nrow(y[!duplicated(y),]) > 1 & nrow(y) >2)) {kmfit = kmeans(y,2)}# Rows are clustered into two centers
    else{kmfit = kmeans(y,1)}
  }
  else{
    if(length(which(!duplicated(y))) > 2| (length(which(!duplicated(y))) > 1 & nrow(y) >2)) {kmfit = kmeans(y,2)}# Rows are clustered into two centers
    else{kmfit = kmeans(y,1)}
  }

  # Find distance between cluster and best value rep(1, number of samples)
  kmcenter = kmfit$centers
  if(nrow(kmcenter) >1){
    drank = sapply(1:nrow(kmcenter), function(i) {
                    distdf = rbind(kmcenter[i,],rep(1,ncol(kmcenter)))
                    dist(distdf)
                    })
    drankdf = data.frame(cluster = c(1:nrow(kmcenter)), drank = drank)
  }
  else{
    drankdf = data.frame(cluster = 1, drank = 1)
  }

  # Find distance of each feature from best performance
  dimp = sapply(1:nrow(y), function(i) {
    distdf = rbind(as.numeric(y[i,]),rep(1,ncol(y)))
    dist(distdf)
  })

  featrank = data.frame(variable = rownames(y), imp = 1/as.numeric(dimp), cluster = kmfit$cluster, stringsAsFactors = F)
  featrank = merge(featrank, drankdf, by = "cluster", all = T)
  featrank$rank = rank(featrank$drank)
  # str(featrank)
  featrank[,c("cluster", "drank")] = NULL
  return(featrank)
}


#'@export
unsuprankagg = function(ranklist, rulemethod = c("cormeth"), perfvar = NA){

  if(!is.na(perfvar)){ranklist = ranklist[ ,!names(ranklist) %in% perfvar]}

  ruletech = c(cormeth = cormeth)

  featres = ruletech[[rulemethod]](ranklist)

  featwt = featres@solution
  featrank = data.frame(variable = names(ranklist), imp = as.numeric(featwt), stringsAsFactors = F)
  featrank$rank = rank(featrank$imp)
  return(featrank)
}

cormeth = function(ranklist, ...){

  # Convert ranklist to actual ranks
  act_rank = apply(ranklist,1,function(x) {x = abs(x);
                                           y = bestNormalize::orderNorm(x[!is.na(x)])$x.t;
                                           newx = x;
                                           newx[!is.na(newx)] = y;
                                           newx})
  act_rank = data.frame(t(act_rank))

  # Optimization: GA
  gafun = function(x){
    # allcor = sapply(1:nrow(act_rank), function(y){cor(x, as.numeric(act_rank[y,]), use = "complete.obs", method = "spearman")
    # })
    allcor = cor(x, t(act_rank), use = "pairwise.complete.obs", method = "kendall")
    fitval = sum(allcor, na.rm = T)
    return(fitval)
  }

  GA = GA::ga(type = "permutation", fitness = gafun, lower = 1, upper = ncol(ranklist), seed = 3, popSize = 100, elitism = 0.1, maxiter = 1000, run = 50, pcrossover = 0.8, pmutation = 0.2, monitor = F)

  return(GA)
}

#'@export
suprankagg = function(ranklist, rulemethod = c("rf", "ridge", "lasso", "rfauto", "autoann"), ytreat = c("No", "rank", "ranknorm", "bellinput"), perfvar = NA){
  if(is.na(perfvar)){print("Cannot perform supervised learning")}

  # Handling Missing Data: replace NA with 0
  ranklist[is.na(ranklist)] = 0

  # str(ranklist)
  # Performnace transformation
  if(ytreat == "rank"){ranklist[,perfvar] = rank(ranklist[,perfvar])}
  else if(ytreat == "ranknorm"){ranklist[,perfvar] = bestNormalize::orderNorm(ranklist[,perfvar])$x.t}
  else if(ytreat == "bellinput"){
    ranklist[,!names(ranklist) %in% perfvar] = apply(ranklist[,!names(ranklist) %in% perfvar], 2,function(x){
      # print(x)
      # if(length(unique(x[!is.na(x)]))>1){z = -1*abs(scale(x))}
      if(length(unique(x[!is.na(x)]))>1){z = -1*abs(x-mean(x[x != 0]))}
      else{z = x*0}
      # print(z)
      z
    })
  }
  # str(ranklist)
  # print(ranklist[is.infinite(ranklist[,perfvar]),perfvar])

  fit = emb_mod(df = ranklist[!is.infinite(ranklist[,perfvar]),], outvar = perfvar, modeltype = "regression", technique = rulemethod)

  featrank = loopdf(df = fit, tech = rulemethod)

  names(featrank) = c("variable", "imp")
  featrank$imp = abs(featrank$imp)
  featrank$rank = rank(1/featrank$imp)

  return(featrank)
}
