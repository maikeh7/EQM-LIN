####################################################################################################
# this code is ONLY for finding the threshold at which the switch from EQM to a linear 
# correction occurs. The threshold is selected based on the mean minimum mean absolute error
# over a K-fold cross validation (see calcMetrics() for how MAE is calculated)
# Please read over the code prior to running--specifically, code under **IMPORTANT** headings below
# should be read thoroughly

##########################################################################################
# Please look at the associated files cvdf_consecYears.csv and RCP_PRCP_1976_2005.Rds" 
# to understand data format
# Additionally, you will need MAE() and calcMetrics() functions in CalcPerformance.R
# cvdf ("cvdf_consecYears.csv ") is a csv file containing variables YEAR and K. YEAR should span the years within 
# the historical period (e.g. 1976-2005). KFOLD denotes the fold for cross validation. Here,
# cvdf is set up for a 5-fold cross-val. Please adjust 'cvdf' according to your needs and 
# preferences. 
# df, RCP_PRCP_1976_2005.Rds" ,is the csv file containing daily precip (mm) data for observations ('PRCP') and raw model data
# downscaled to station locations ('test.pred.mu') 
# other variables in 'df' that you may need to change based on your data:
# Station_ID: ID for climate station
# MONTH: month of year (1-12)
# YEAR: 4 digit integer value for year
# WRF_PRCP: raw WRF. This is the nearest WRF grid point to a station location
# LAT, LON , LON_LCC, LAT_LCC are lat/lon
# PRCP: station (observed precip) (mm)
# test.pred.mu: WRF precip downscaled to station locations (mm)
# by Maike Holthuijzen
# updated 6/27/2022
##################################################################################################

library(dplyr)
library(qmap)
library(reshape2)
library(data.table)

# example data
df = readRDS("RCP_PRCP_1976_2005.Rds")
cvdf = read.csv("cvdf_ConsecYears.csv")

# grid search values for threshold
Qthreshes = seq(.60, .95, by = .01)

maeList = list()
# vector to hold MAE95 results
meanMAEVec = vector(length = length(Qthreshes))

# vector to hole MAE results 
MAEvec = vector(length = length(Qthreshes))
for (m in 1:length(Qthreshes)){
  #print(m)
  Qthresh = Qthreshes[m]
  CVresDF = data.frame()
  kFolds = 5
  for (i in 1:kFolds) {
    print(i)
    trainYears = filter(cvdf, K != i)$YEAR
    testYears = filter(cvdf, K == i)$YEAR
    traindf = filter(df, YEAR %in% trainYears)
    
    # for each fold, the threshold is calculated. This determines the split between
    # EQM and the linear correction
    myquant = quantile(traindf$test.pred.mu[traindf$test.pred.mu > 0], probs = Qthresh)
    testdf = filter(df, YEAR %in% testYears)
    LinCorrDF = data.frame()
    # loop thru each month and perform EQM-LIN based on the current value of Qthresh
    for (j in 1:12){
      temp = filter(traindf, MONTH == j)
      
      obsall = temp$PRCP
      modall = temp$test.pred.mu
      
      temptest = filter(testdf, MONTH == j)
      obstestall = temptest$PRCP
      modtestall = temptest$test.pred.mu
      
      # this is the portion of the test set to be corrected via EQM
      modtesteqmDF = filter(temptest, test.pred.mu < myquant)
      modtesteqm = modtesteqmDF$test.pred.mu
      
      # this is the portion of the test set to be corrected by linear correction
      modtestLinearDF = filter(temptest, test.pred.mu >= myquant)
      modtestLinear = modtestLinearDF$test.pred.mu
      
      #IMPORTANT READ ####################################################################################
      # construct EQM TF (use all data here)
      # wet.day should be adjusted as desired! 
      # **NOTE: using wet.day = F will result in different results
      # if there is no lower threshold, use wet.day = FALSE (see documentation for 
      # fitQmapQUANT in qmap package)**
      # In paper, I used wet.day = F and then implemented threshold, but it is likely better
      # to use a wet day adjustment within quantile mapping if your precip data has a lower precip threshold
      # NOTE--if you use wet.day = TRUE, I found that for Qthresh < .6, MAE95 continues to 
      # decrease until 0.35 **BUT** resulted in an increase in MAE
      # Therefore, I would recommend the lower bound on Qthresh should be > 0.6 (e.g Qthreshes = seq(.60, .95, by = .01) )
      # #####################################################################################################
      
      # this is what is used in Holthuijzen et al 2022. 
      tf = fitQmapQUANT(obsall, modall, qstep = 0.0001, wet.day = FALSE)
      
      # IMPORTANT:
      # Slightly better results MAY be obtained using wet.day = 0.1) (my precip data has a 0.1 mm min threshold)
      # comment out whichever option you do not use
       tf = fitQmapQUANT(obsall, modall, qstep = 0.0001, wet.day = 0.1)
      
      
      mydat = data.frame(modq = tf$par$modq, fitq = tf$par$fitq)
      
      # mydathi is the portion of the EQM TF containing model data greater than 
      # or equal to myquant
      mydathi = filter(mydat, modq >= myquant)[1, ]
      
      # constantCorr determines the constant linear correction (difference between
      # observed and model data right where the threshold is)
      constantCorr = mydathi$fitq - mydathi$modq
      
      # apply EQM to model data in 'EQM' test set  (lower quantiles)
      eqmcorr = doQmapQUANT(modtesteqm, tf, type = "tricub")
      modtesteqmDF$LinearCorr = eqmcorr
      
      # apply linear correction (upper quantiles of model data plus constantCorr)
      linearCorr = modtestLinear + constantCorr
      modtestLinearDF$LinearCorr = linearCorr
     
      # some plotting code if you want to plot
      # mydathidf = filter(mydat, modq >= myquant) 
      # mydatlodf = filter(mydat, modq < myquant)
      # mydatlodf$ynew = mydatlodf$fitq
      
      # mydathidf$ynew = mydathidf$modq + constantCorr
      # newdf = rbind(mydatlodf, mydathidf)
      # newdf$Quantile = seq(0,1,by=0.0001)
      # newdf$f1 = newdf$fitq
      # newdf$f2 = newdf$modq + constantCorr
      # plot(newdf$modq, newdf$ynew, pch = 19, main = paste("kfold ", i, "MONTH ", j), cex=.5, xlim = c(0,40))
      # abline(v=myquant, col = "blue")
      
      # bind the EQM- and linear-corrected data into one dataframe
      AllcorrDF = rbind(modtesteqmDF, modtestLinearDF)
      
      # add to the big dataframe containing data for all K folds 
      LinCorrDF = rbind(AllcorrDF, LinCorrDF)
    }
    
    # these variable names will need to be adjusted for your data! See above for variable definitions
    LinCorrDF = dplyr::select(LinCorrDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, LinearCorr)
    LinCorrDF = arrange(LinCorrDF ,YEAR ,MONTH, DAY , Station_ID)
    kfoldDF = LinCorrDF
    varname = "EQMLIN"
    kfoldDF[[varname]] = kfoldDF$LinearCorr
    
    
    #kfoldDF = dplyr::select(kfoldDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, all_of(varname))
    kfoldDF = dplyr::select(kfoldDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, all_of(varname))
    
    kfoldDF$KFOLD = i
    #CVresDF is the dataframe containing results for all folds
    CVresDF = rbind(kfoldDF, CVresDF)
  }
  
  # Now, calculate performance metrics
  # Need calcMetrics() function for this!
  metrics1 = calcMetrics(CVresDF)
  mylist = metrics1$mylist
  
  biglist = rbindlist(mylist) %>% dplyr::select("maeRaw", "maeCorr", "mae95Raw", "mae95Corr","MONTH", "KFOLD")
  #biglist = rbindlist(mylist) %>% dplyr::select("maxMonthlyMod", "maxMonthlyObs", "maxMonthlyCorr",  "MONTH", "KFOLD")
  meltdf = melt(biglist, id.vars = c("MONTH", "KFOLD"))
  
  # 'variable' refers to the performance metric (e.g. maeRaw, maeCorr, etc)
  # 'value' is the value of the performance metric
  # here, we are interested in minimizing MAE95--> 'mae95Corr'
  res2 = meltdf %>% group_by(variable) %>% summarise(meanMAE = mean(value))
  maeList[[m]] = res2
  meanMAEVec[m] = res2$meanMAE[4] #this is MAE95 for EQM-LIN corrected data (**the only metric we are interested in here**)
  MAEvec[m] = res2$meanMAE[2] # this is MAE for EQM-LIN corrected data (can inspect this as well)
  }

# plot the vector containing mean MAE95 results against the possible values for thresholds
# choose the value of Qthresh that results in minimum MAE95
plot(Qthreshes, meanMAEVec)
Qthreshes(which.min(meanMAEVec))
