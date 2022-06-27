
library(data.table) # for data.table::rbindlist
library(reshape2) # for reshape2::melt
library(dplyr) # for data manipulation
library(ggplot2)


# function to calculate mean absolute error (MAE)
MAE = function(obs, mod){
  mean(abs(obs - mod))
}

###############################################################################
# function to calculate performance metric for the data object resulting from 
# cross-validation in EQMLIN.R
# See example below
# dat is the resulting dataframe and should look as follows:
# nfolds are the number of folds in the cross-validation. Unless modified, it should be 5
calcMetrics = function(dat, nfolds = 5){
  nfolds = length(unique(dat$KFOLD))
  mylist=list()
  # VAR refers to the bias-correction method being used. I used this code to test other bias correction methods, 
  # but for EQM-LIN, it will be 'EQM_LIN'
  VAR = colnames(dat)[7]
  mydf = data.frame()
  
  # for each of the 5 folds, calculate performance metrics
  for (k in 1:nfolds){
    myfold = k
    kdf = filter(dat, KFOLD == myfold)
    
    # get estimated quantiles
    # obs = observed data, mod = raw model data, corrP = bias-corrected model data
    mod = kdf$test.pred.mu
    qmod = quantile(mod, probs = seq(0,1, by = 0.0001))
    obs = kdf$PRCP
    qobs = quantile(obs, probs = seq(0,1, by = 0.0001))
    corrP = kdf[[VAR]]
    qcorrP = quantile(corrP, probs = seq(0,1,  by = 0.0001))
    
    # we only want values greater than 0 mm precip
    corrP = corrP[corrP > 0]
    obs = obs[obs > 0]
    
    #calculated for wet day precip! Modify as necessary 
    mod = ifelse(mod < .1, 0 , mod)
    mod = mod[mod > 0]
    
    # get 95th quantile for MAE95
    q95o = quantile(obs, probs = .95)
    q95w = quantile(mod, probs = .95)
    q95c = quantile(corrP, probs = .95)
    
    # note that the values above the 95th quantile for all 3 data types may not be the same. 
    # this is taken care of below
    mymod = mod[mod >= q95w]
    myobs = obs[obs >= q95o]
    mycorr = corrP[corrP >= q95c]

    # get maximum length of values above 95th quantile for each of the 3 data types
    # Then, estimate the extreme tail quantiles. Number of extreme tail quantiles will be equal to 'mymax'
    mymax = max(length(mymod), length(myobs), length(mycorr))
    q95corrP = quantile(corrP, probs = seq(.95, 1, length.out = mymax))
    q95obs = quantile(mod, probs = seq(.95, 1, length.out = mymax))
    q95mod = quantile(obs, probs = seq(.95, 1, length.out = mymax))
    
    # calculate MAE and MAE 95
    maeRaw = MAE(qobs, qmod)
    maeCorr = MAE(qobs, qcorrP)
    mae95Raw = MAE(q95obs, q95mod)
    mae95Corr = MAE(q95obs, q95corrP)
    
    # This is how I named the metrics. Adjust if you want.
    # Raw = raw model data, Corr = corrected model data. 'myfold' = the fold number
    metrics = c(maeRaw, maeCorr, mae95Raw, mae95Corr, myfold)
    
    # Overall metrics (not by month)
    mydf = rbind(mydf, metrics)
    
    # calculate metrics by month--this is what is used in the paper
    kfoldDF = data.frame()
    for (j in 1:12){
      monthdf = filter(kdf, MONTH == j)
      mod = monthdf$test.pred.mu
      qmod = quantile(mod, probs = seq(0,1, by = 0.0001))
      obs = monthdf$PRCP
      qobs = quantile(obs, probs = seq(0,1, by = 0.0001))
      corrP = monthdf[[VAR]]
      qcorrP = quantile(corrP, probs = seq(0,1, by = 0.0001))
      
      
      
      
      corrP = corrP[corrP > 0]
      obs = obs[obs > 0]
      mod = ifelse(mod < .1, 0 , mod) #calculated for wet day precip
      mod = mod[mod > 0]
      
      q95o = quantile(obs, probs = .95)
      q95w = quantile(mod, probs = .95)
      q95c = quantile(corrP, probs = .95)
      
      mymod = mod[mod >= q95w]
      myobs = obs[obs >= q95o]
      mycorr = corrP[corrP >= q95c]
      
      mymax = max(length(mymod), length(myobs), length(mycorr))
      q95corrP = quantile(corrP, probs = seq(.95, 1, length.out = mymax))
      q95obs = quantile(mod, probs = seq(.95, 1, length.out = mymax))
      q95mod = quantile(obs, probs = seq(.95, 1, length.out = mymax))
      
      # I do not use these in the paper. Was more exploratory
      maxMonthlyObs = max(obs)
      maxMonthlyMod = max(mod)
      maxMonthlyCorr = max(corrP)
      
      maeRaw = MAE(qobs, qmod)
      maeCorr = MAE(qobs, qcorrP)
      mae95Raw = MAE(q95obs, q95mod)
      mae95Corr = MAE(q95obs, q95corrP)
      metrics = c(maeRaw, maeCorr, mae95Raw, mae95Corr, maxMonthlyMod, maxMonthlyObs, maxMonthlyCorr, j, myfold)
      
      
      kfoldDF = rbind(kfoldDF, metrics)
      
    }
    # max monthly metrics are not used in paper. They were for exploratory analyses
    names(kfoldDF) = c("maeRaw", "maeCorr", "mae95Raw", "mae95Corr", "maxMonthlyMod", 
                       "maxMonthlyObs", "maxMonthlyCorr",  "MONTH", "KFOLD")
    mylist[[k]] = kfoldDF
  }
  names(mydf) = c("maeRaw", "maeCorr", "mae95Raw", "mae95Corr", "KFOLD")
  return(list(mylist = mylist, mydf = mydf))
}

#############################################
# EXAMPLE 
# read in resulting dataframe from EQMLIN.R
#############################################
# this is the result from using wet.day = FALSE in FindThreshold.R (used in paper)
# the threshold is 0.79 
EQMLIN = readRDS("EQMLIN79.Rds")

# You can compare to this, which is the result when wet.day = TRUE in FindThreshold.R
# the threshold is 0.60 here
# EQMLIN = readRDS("EQMLIN6.Rds")

# calculate MAE/MAE95 
# note that results are store PER FOLD. For cross-validated metrics, average over folds
resLIN = calcMetrics(EQMLIN, nfolds = 5)

# These are the results by month
mylist = resLIN$mylist

# Some quick manipulation to get a dataframe containing month and value of metrics
biglist = rbindlist(mylist) %>% dplyr::select("maeRaw", "maeCorr", "mae95Raw", "mae95Corr","MONTH", "KFOLD")

meltdf = melt(biglist, id.vars = c("MONTH", "KFOLD"))



