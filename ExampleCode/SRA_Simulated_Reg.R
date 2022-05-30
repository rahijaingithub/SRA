# Ensemble Rank Aggregation: Regression
plist = c(75)
nlist = c(100)
targetfeat = c(10)
correlationval = c(5)
correlation_var= c(15)
var_effect1 = c(0.9)
var_effect2 = c(0.9)
var_effect3 = c(0.9)
var_effect4 = c(0.9)
bootval = c(300)

paralist = expand.grid(p = plist, n = nlist, mar = targetfeat, cor = correlationval, corvar = correlation_var, ve1 = var_effect1, ve2 = var_effect2, ve3 = var_effect3, ve4 = var_effect4, bootval = bootval)
testcases = paralist

testcases[2,] =c(100,100,10,5,15,0.5,0.5,0.5,0.5, 100)
testcases[3,] =c(175,275,15,5,15,0.4,-0.8,0.4,-0.8, 100)
testcases[4,] =c(75,275,15,5,15,0.4,-0.8,0.4,-0.8, 100)
testcases[5,] =c(75,225,15,5,15,0.4,-0.8,0.4,-0.8, 200)
testcases[6,] =c(125,225,20,5,15,0.4,-0.8,0.4,-0.8, 200)
testcases[7,] =c(100,95,10,5,15,0.5,0.5,0.5,0.5, 100)

techniques = c("ridge")
for(techie in techniques){#nrow(testcases)
  val = pbapply::pblapply(2:2, function(scene){
    cat(scene, " ")
    p = testcases$p[scene]
    n = testcases$n[scene]
    mar = testcases$mar[scene]
    corval = testcases$cor[scene]
    corvar = testcases$corvar[scene]
    var_effect = as.numeric(testcases[scene, c("ve1", "ve2", "ve3", "ve4")])
    bootval = testcases$bootval[scene]
    modeltype = "regression"
    vareffect = stringi::stri_replace_all_fixed(stringi::stri_replace_all_fixed(paste(var_effect, collapse = ""), "0.", ""), "-","neg")
    fullperfres = pbapply::pblapply(1:10, function(j){
      set.seed(j)
      sampleregdf= regdataset(p, train_sample = n, var = "No_Int", setting = "Correlation", datadist = c("normal"), corr = c("fixed"), seed = j, main_var=mar, var_effect=var_effect, correlation_var=corvar, correlation_val=corval)
      traindf = sampleregdf[[1]]
      testdf = sampleregdf[[2]]

      # Run each ensemble method
      tech = techie
      res = ensem_homo(df = traindf, outvar = "y", modeltype = modeltype, technique = tech, sampletype = "bootsample", fstype = "randomfeat", maxfeatsample = ncol(traindf)-1, boots = bootval, savememory = F, Trace = F, parallel = T, smart = T, perfparallel = F)

      alltech = alltechlist(rulemethod = c("sd","min", "max", "mean", "median","freq", "ridge", "rfauto", "lasso", "rra",  "coefv", "ttest", "wmw"))

      perfres = lapply(1:nrow(alltech), function(i){
        res[[2]]$perf = res[[2]]$perf/max(res[[2]]$perf)
        # Rank Aggregation Method
        rameth = c(single = rankaggfun)

        # Rank Aggregation with km based cutoff
        rulemethod = alltech[i,1]
        i = c(rulemethod, tech)
        # print(i)
        if(!rulemethod %in% c("serial", "parallel", "hybrid", "wt_hybrid")){

          if(rulemethod == "cormeth"){
            rankmethod = "unsuprankagg"
            finrank = rankaggfun(ranklist = res[[2]], rankmethod = rankmethod, rulemethod = rulemethod, perfvar = "perf")
          }
          else if(rulemethod %in% c("ridge","rf", "lasso","ann", "rfauto","autoann")){
            rankmethod = "suprankagg"
            ytreat = "No"
            i = c(rulemethod, ytreat)
            finrank = rankaggfun(ranklist = res[[2]], rankmethod = rankmethod, rulemethod = rulemethod, ytreat=ytreat, perfvar = "perf")
          }
          else{
            rankmethod = "rulebase"
            finrank = rankaggfun(ranklist = res[[2]], rankmethod = rankmethod, rulemethod = rulemethod, perfvar = "perf")
          }

          # KM based Cut-off
          kmdf = finrank
          kmdf$rank = NULL
          kmdf = data.frame(t(kmdf))
          names(kmdf) = kmdf[1,]
          kmdf = data.frame((kmdf[-1,]))
          kmdf = as.data.frame(lapply(kmdf, as.numeric))
          kmres = km(y = kmdf)
          finrank = kmres
        }
        else if(rulemethod %in% c("wt_hybrid")){
          finrank = rameth[[rulemethod]](ranklist = res[[2]], rankmethod = NA, rulemethod = c("min", "max", "mean", "median","freq", "ridge","rfauto", "lasso", "rra", "coefv", "ttest", "wmw"), perfvar = "perf", metric = c("featimp"), traindf = traindf, outvar = "y", modeltype = modeltype, wtfactor = 10)
        }
        else{
          finrank = rameth[[rulemethod]](ranklist = res[[2]], rankmethod = NA, rulemethod = c("min", "max", "mean", "median","freq", "ridge","rfauto", "lasso", "rra", "coefv", "ttest", "wmw"), perfvar = "perf", metric = c("featimp"))
        }

        selfeat = outdat(tech = "km", resdf = finrank, pretech = paste("randomfeat", i[1], i[2], sep = "_"), algo = "Ensemble", martarget = mar, perfmetric = "F1")

        selfeat[1,c("perf", "type")]

        # Predictive Performance
        feat = finrank$variable[finrank$rank == min(finrank$rank) & finrank$imp != 0]
        mp = predperf(traindf = traindf, testdf = testdf, outvar = "y", modeltype = modeltype, selfeat = feat, pretech = paste(i[1], i[2], sep = "_"), martarget = mar, interaction = F)
        # print(mp)
        out = c(selfeat[1,c("perf", "type")],mp[, c("pp", "target", "noise")])
        out = data.frame(out)
      })
      perfres
    })
    fullperfres2 = purrr::flatten(fullperfres)
    fullperfdf = do.call(rbind, fullperfres2)

    filename = paste("10trial_reg_ensemble_pn_",p,n,"_mainvar_",mar,"_vareffect_",vareffect, "_corvar_", corvar, "_corval_", corval, "_",techie,"__rankagghyb_",bootval,"_randomfeat_kmcutoff.csv", sep = "")
    write.csv(fullperfdf, filename)

    summarydf = aggregate(. ~ type, data = fullperfdf, FUN = function(x) c(mean = mean(x), ci = unlist(colcint(x))) )
    filename = paste("reg_ensemble_pn_",p,n,"_mainvar_",mar,"_vareffect_",vareffect, "_corvar_", corvar, "_corval_", corval, "_",techie,"__rankagghyb_",bootval,"_randomfeat_kmcutoff.csv", sep = "")
    write.csv(summarydf, filename)
  })
}
