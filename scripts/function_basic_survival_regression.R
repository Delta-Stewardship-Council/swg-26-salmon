## This function performs a basic regression for our initial survival analysis,
##   of the form y = f(x1, x2, ...), where y is survival/passage and x are environmental/etc. variables
## Inputs:
##     data: R dataframe, where each row should represent a single fish
##        y: name of the column representing survival/passage data (1 if fish survived, 0 otherwise)
##        x: vector, name of column(s) for dependent variables to include in the model (flow, temp, etc.)
## mod_type: type of regression model. Currently implemented options include "LR" and "GAM"
##          "GLM": classic binary logistic regression, using base R
##          "GAM": generalized additive model, using mgcv package
## Outputs:
##       res: raw model output of glm/gam
## preds_out: predictions of model for plotting

require(mgcv)

basic_survival_regression = function(data, y, x, mod_type = "GLM", resp_spec = NULL){
  
  ## make new data frame for predictions
  if(is.null(resp_spec)){
    nunique = function(k){length(unique(k))}
    resp_spec = ifelse(sapply(data[, x], nunique) > 5 & sapply(data[, x], is.numeric), "cont", "fac")
  }
  
  xvals = list()
  for(i in 1:length(x)){
    xvar = x[i]
    if(resp_spec[i] == "cont"){
      xvals[[xvar]] = seq(min(data[, xvar]), max(data[, xvar]), length.out = 25)
    }else if(resp_spec[i] == 'fac'){
      xvals[[xvar]] = unique(data[, xvar])
    }
  }
  
  newdata = do.call(expand.grid, xvals)
  
  ## fit models
  if(mod_type == "GLM"){
    #get formula
    glm_form = as.formula(paste(y, "~", paste(x, collapse = " + ")))
    
    #fit model and make predictions
    res = glm(glm_form, data = data, family = binomial(link = 'logit'))
    preds = predict(res, newdata = newdata, type = "response", se.fit = T)
    preds_out = cbind(prediction = preds$fit, se = preds$se.fit, newdata)
  }else if(mod_type == "GAM"){
    
    #get formula
    gam_form = paste(y, "~ ")
    for(i in 1:length(x)){
      xvar = x[i]
      if(i > 1){ gam_form = paste0(gam_form, " + ")}
      if(resp_spec[i] == "cont"){
        gam_form = paste0(gam_form, "s(", xvar, ")")
      }else if(resp_spec[i] == 'fac'){
        gam_form = paste0(gam_form,  xvar)
      }
    }
    gam_form = as.formula(gam_form)
    #gam_form = as.formula(paste(y, "~", paste(paste0("s(", x, ")"), collapse = " + ")))
    
    #fit model and make predictions
    res = gam(gam_form, 
              data = data, 
              family = binomial(link = 'logit'))
    preds = predict(res, newdata = newdata, type = "response", se.fit = T)
    preds_out = cbind(prediction = preds$fit, se = preds$se.fit, newdata)
  }
  
  return(list(res = res, predictions = preds_out))
  
}
