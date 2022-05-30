#' This function runs embedded model
#'
#'@export
emb_mod = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", ...){
 alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
             step = NA, ann = anntech, dt = NA, reg= NA, rfauto = autorftech, flasso = flasso, falasso = falasso, fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

 if(length(list(...)) != 0){
         argvec = c(...)
         if (technique %in% c("step", "dt")){fit = NA}
         else{
                 if(grepl("auto", technique)){technique = stringi::stri_replace_all_fixed(technique, "auto","")}
                 fit = do.call(alltech[[technique]],c(argvec, list(df = df, outvar=outvar, modeltype = modeltype)))}
 }
 else{
         if (technique %in% c("step", "dt")){fit = NA}
         else{fit = alltech[[technique]](df = df, outvar = outvar, modeltype = modeltype)}
 }
 return(fit)
}

#'@export
emb_mod_int = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", ...){
        alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
                    step = NA, ann = anntech, dt = NA, reg= NA, rfauto = autorftech, flasso = flasso, falasso = falasso, fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

        if(length(list(...)) != 0){
                argvec = c(...)
                if (technique %in% c("step", "dt")){fit = NA}
                else{
                        if(grepl("auto", technique)){technique = stringi::stri_replace_all_fixed(technique, "auto","")}
                        fit = do.call(alltech[[technique]],c(argvec, list(df = df, outvar=outvar, modeltype = modeltype)))}
        }
        else{
                if (technique %in% c("step", "dt")){fit = NA}
                else{fit = alltech[[technique]](df, outvar, modeltype, interaction = T)}
        }
        return(fit)
}


#'@export
emb_mod_cov = function(df = traindf, outvar = "y", modeltype = "regression", technique = "rf", covs= NA, ...){
        alltech = c(rf = rftech, lasso = lasso, alasso = alasso, ridge = ridge, aridge = aridge, enet = enet, aenet = aenet,
                    step = NA, ann = anntech, dt = NA, reg= NA, rfauto = autorftech, flasso = flasso, falasso = falasso, fridge = fridge, faridge = faridge, fenet = fenet, faenet = faenet, autoann = autoanntech)

        if(length(list(...)) != 0){
                argvec = c(...)
                if (technique %in% c("step", "dt")){fit = NA}
                else{
                        if(grepl("auto", technique)){technique = stringi::stri_replace_all_fixed(technique, "auto","")}
                        fit = do.call(alltech[[technique]],c(argvec, list(df = df, outvar=outvar, modeltype = modeltype)))}
        }
        else{
                if (technique %in% c("step", "dt")){fit = NA}
                else{fit = alltech[[technique]](df, outvar, modeltype, covs = covs)}
        }
        return(fit)
}
