#########################################################################################3
# This code implements EQM-LIN in a 5-fold cross-validation setting.
# You should FIRST RUN FindThreshold.R to get the threshold!!

# The code: loop over each of the 5 folds and month and perform EQM-LIN. Results for all folds 
# are stored in CVresDF, a data.frame
# NOTE: you may wish to adjust qstep or the value for wet.day in fitQmapQUANT()!
# To make it easier for others to modify, I did not create a function for the implementation
# of EQM-LIN below. Please modify as necessary for your work.
# You will have to modify this code so that it works for your data/variable names 
# Additionally, you will need MAE() and calcMetrics() functions in CalcPerformance.R
##########################################################################################
# Please look at the associated files cvdf_consecYears.csv and RCP_PRCP_1976_2005.Rds" 
# to understand data format
# cvdf ("cvdf_consecYears.csv ") is a csv file containing variables YEAR and KFOLD. YEAR should span the years within 
# the historical period (e.g. 1976-2005). KFOLD denotes the fold for cross validation. Here,
# cvdf is set up for a 5-fold cross-val. Please adjust 'cvdf' according to your needs and 
# preferences. 
# df, RCP_PRCP_1976_2005.Rds" ,is the csv file containing daily precip (mm) data for observations ('PRCP') and raw model data
# downscaled to station locations ('test.pred.mu') 
# other variables in 'df' that you may need to change based on your data:
# Station_ID: ID for climate station
# MONTH: month of year (1-12)
# YEAR: 4 digit integer value for year
# by Maike Holthuijzen
# updated 6/27/2022
##########################################################################################

cvdf = read.csv("cvdf_consecYears.csv")
df= readRDS("Data/RCP_PRCP_1976_2005.Rds")

CVresDF = data.frame()
Qthresh = .79 # This needs to be adjusted based on results from grid search (FindThreshold.R)!
for (i in 1:kFolds) {
  print(i)
  trainYears = filter(cvdf, K != i)$YEAR
  testYears = filter(cvdf, K == i)$YEAR
  traindf = filter(df, YEAR %in% trainYears)
  myquant = quantile(traindf$test.pred.mu[traindf$test.pred.mu > 0], probs = Qthresh)
  testdf = filter(df, YEAR %in% testYears)
  LinCorrDF = data.frame()
  for (j in 1:12){
    temp = filter(traindf, MONTH == j)
    obsall = temp$PRCP
    modall = temp$test.pred.mu
    
    temptest = filter(testdf, MONTH == j)
    obstestall = temptest$PRCP
    modtestall = temptest$test.pred.mu
    #qqtest = fitQmapQUANT(obstestall, modtestall, qstep = 0.0001, wet.day = 0.1)
    #plot(qqtest$par$modq, qqtest$par$fitq, xlab = "WRF quantiles", ylab = "Obs quantiles", main = "QQ map Jan (Test set)", pch = 19)
    #abline(0,1, col = "blue")
    #myquant = quantile(modall[modall > 0], probs = Qthresh)
    
    modtesteqmDF = filter(temptest, test.pred.mu < myquant)
    modtesteqm = modtesteqmDF$test.pred.mu
    
    modtestLinearDF = filter(temptest, test.pred.mu >= myquant)
    modtestLinear = modtestLinearDF$test.pred.mu
    
    # important to keep wet.day == TRUE, adjust as necessary
    # If your precip data does not have a lower threshold, use wet.day = FALSE
    tf = fitQmapQUANT(obsall, modall, qstep = 0.0001, wet.day = 0.1)
    mydat = data.frame(modq = tf$par$modq, fitq = tf$par$fitq)
    
    mydathi = filter(mydat, modq >= myquant)[1, ]
    mydathidf = filter(mydat, modq >= myquant)
    mydatlodf = filter(mydat, modq < myquant)
    
    constantCorr = mydathi$fitq - mydathi$modq
    
    eqmcorr = doQmapQUANT(modtesteqm, tf, type = "tricub")
    modtesteqmDF$LinearCorr = eqmcorr
    
    linearCorr = modtestLinear + constantCorr
    modtestLinearDF$LinearCorr = linearCorr
    # bind the EQM- and linear-corrected data into one dataframe
    AllcorrDF = rbind(modtesteqmDF, modtestLinearDF)
    
    # add to the big dataframe containing data for all K folds 
    LinCorrDF = rbind(AllcorrDF, LinCorrDF)
  }
  # this will need to be adjusted for your data! See above for variable definitions
  LinCorrDF = dplyr::select(LinCorrDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, LinearCorr)
  LinCorrDF = arrange(LinCorrDF ,YEAR ,MONTH, DAY , Station_ID)
  kfoldDF = LinCorrDF
  varname = "EQMLIN"
  kfoldDF[[varname]] = kfoldDF$LinearCorr
  
  
  #kfoldDF = dplyr::select(kfoldDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, all_of(varname))
  kfoldDF = dplyr::select(kfoldDF, YEAR, MONTH, DAY, Station_ID, PRCP, test.pred.mu, varname)
  
  kfoldDF$KFOLD = i
  #CVresDF is the dataframe containing results for all folds
  CVresDF = rbind(kfoldDF, CVresDF)
}

CVresDF = readRDS("CVresLinDF.Rds")
metrics1 = calcMetrics(CVresDF)
mylist = metrics1$mylist

# see below for definitions of metrics
biglist = rbindlist(mylist) %>% dplyr::select("maeRaw", "maeCorr", "mae95Raw", "mae95Corr","MONTH", "KFOLD")

#biglist = rbindlist(mylist) %>% dplyr::select("maxMonthlyMod", "maxMonthlyObs", "maxMonthlyCorr",  "MONTH", "KFOLD")

meltdf = melt(biglist, id.vars = c("MONTH", "KFOLD"))
#######################################################3
#"maeRaw", "maeCorr", "mae95Raw", "mae95Corr"
# maeRaw = MAE for raw model data
# maeCorr = MAE for EQM-LIN corrected model data
# mae95Raw = MAE for raw model data
# mae95Corr = MAE for EQM-LIN corrected model data