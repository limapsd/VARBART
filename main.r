#################################################################################################
#
#  MINESSOTA BART - PROTOTYPE
#
#################################################################################################
# Author : Pedro Antonio Sa Barreto de Lima
# The University of Texas at Austin - McCombs Business School
# Email : pedro.lima@mccombs.utexas.edu
#################################################################################################
# This code is trying to be based on the Carriere et al(2019) Matlab implementation

# Loading required libraries:
library("MASS")
library("stochvol")
library("dbarts")
library("bayesianVARs")

source("varbartfunc.R")
source("flexBART1.R")

################# Data and Manipulations ################:

variables <- c("GDPC1", "CPIAUCSL", "INDPRO", "S&P 500")
#variables  = colnames(usmacro_growth)[1:13]   

train_data     = 100 * usmacro_growth[1:231, variables]
test_data      = 100 * usmacro_growth[232:234, variables]
lags           = 2
h              = 3


############# Models #############:

index = which(variables %in% "CPIAUCSL")

#################### My Implementation #####################

# BART without SV:
bart.h = VARBART(train_data, p =lags, fhorz = h, sv = FALSE)
get_barth.forecast = bart.h$fcst
bart.h.forecast = get_barth.forecast[,,"CPIAUCSL"]

quantile.bart.h = apply(bart.h.forecast, 2 , quantile, probs = c(0.05,0.5,0.95))
point.wise.bart.h = quantile.bart.h["50%",]

# BART with SV
bart.sv = VARBART(train_data, p =lags, fhorz = h, sv = TRUE)
get_bartsv.forecasts = bart.sv$fcst
bart.forecast.sv = get_bartsv.forecasts[,,"CPIAUCSL"]


quantile.bart.sv = apply(bart.forecast.sv, 2 , quantile, probs = c(0.05,0.5,0.95))
point.wise.bart.sv = quantile.bart.sv["50%",]
###########################################################

########### Observed and Fitted ################### 
index = "CPIAUCSL"

Y = bart.h$Y
ylimb = c(-2, 3)
par(mfrow = c(1,2))

plot(Y[,index], col=1 , type = "l", ylab = colnames(Y)[index], xlab = "Time", main = "BART")
lines(bart.h$Y.fit.BART[,index] , col = 2)

plot(Y[,index], col=1 , type = "l", ylab = colnames(Y)[index], xlab = "Time", main = "BART-SV")
lines(bart.sv$Y.fit.BART[,index] , col = 2)

############ Predictions ###########################
plot(test_data[1:h,index], col=1 , type = "l", ylab = index, xlab = "Time", main = "BART", ylim = ylimb)
lines(quantile.bart.h["50%",] , col = 2)
lines(quantile.bart.h["5%",], lty = "dashed", col =2)
lines(quantile.bart.h["95%",], lty = "dashed", col =2)

plot(test_data[1:h,index], col =1, type = "l", ylab = index, xlab = "Time",main = "BART-SV", ylim = ylimb)
lines(quantile.bart.sv["50%",] , col =2)
lines(quantile.bart.sv["5%",], lty = "dashed", col =2)
lines(quantile.bart.sv["95%",], lty = "dashed", col =2)


########### Graphs and Analyzes #################

Forecast.Comparison.Tbl = matrix(NA, nrow = 5, ncol = fhz)
row.names(Forecast.Comparison.Tbl) = c("BART","BART.SV","BVAR","BVAR-MINN","BVAR-MINN-SV")

errors.m = matrix(NA, fhz, ncol =5 )
rownames(errors.m) = names(point.wise.bart.h)
colnames(errors.m) = c("BART","BART.SV","BVAR","BVAR-MINN","BVAR-MINN-SV")

errors.m[,1] = (point.wise.bart.h - test_data[1:fhz,var.int])^2
errors.m[,2] = (point.wise.bart.sv - test_data[1:fhz,var.int])^2
errors.m[,3] = (forecast.1 - test_data[1:fhz,var.int])^2
errors.m[,4] = (forecast.3 - test_data[1:fhz,var.int])^2
errors.m[,5] =  (forecast.2 - test_data[1:fhz,var.int])^2

mse.errors.m = matrix(NA, fhz, ncol =5)
rownames(mse.errors.m) = names(point.wise.bart.h)
colnames(mse.errors.m) =  c("BART","BART.SV","BVAR","BVAR-MINN","BVAR-MINN-SV")

for(i in 1:5){
  for(j in 1:fhz){
    if(j == 1){
      print(errors.m[j,i])
      mse.errors.m[j,i] = errors.m[j,i]
    }else{
      mse.errors.m[j,i] = mean(errors.m[1:j,i])
    }
  }
}



########## Forecasting Table ###############

Forecast.Comparison.Tbl[1,] = point.wise.bart.h
Forecast.Comparison.Tbl[2,] = point.wise.bart.sv
Forecast.Comparison.Tbl[3,] = forecast.1
Forecast.Comparison.Tbl[4,] = forecast.3
Forecast.Comparison.Tbl[5,] = forecast.2

t(Forecast.Comparison.Tbl)/(Forecast.Comparison.Tbl[2,])

######## Forecasting plots ###############


par(mfrow = c(2,2))

ylimb = c(-2, 3)


plot(test_data[1:fhz,var.int], col=1 , type = "l", ylab = var.int, xlab = "Time", main = "BART", ylim = ylimb)
lines(point.wise.bart.h , col = 2)
lines(quantile.bart.h["5%",], lty = "dashed", col =2)
lines(quantile.bart.h["95%",], lty = "dashed", col =2)

plot(test_data[1:fhz,var.int], col =1, type = "l", ylab = var.int, xlab = "Time",main = "BART-SV", ylim = ylimb)
lines(point.wise.bart.sv , col =2)
lines(quantile.bart.sv["5%",], lty = "dashed", col =2)
lines(quantile.bart.sv["95%",], lty = "dashed", col =2)

plot(test_data[1:fhz,var.int], col=1 , type = "l", ylab = var.int, xlab = "Time", main = "BVAR", ylim = ylimb )
lines(forecast.1, col = 2)
lines( summary(pred.1)$prediction_quantiles[1,,index], lty = "dashed", col =2)
lines( summary(pred.1)$prediction_quantiles[3,,index], lty = "dashed", col =2)

plot(test_data[1:fhz,var.int], col =1, type = "l", ylab = colnames(test_data)[2], xlab = "Time",main = "BVAR-Minn-SV", ylim = ylimb)
lines(forecast.2, col =2)
lines( summary(pred.2)$prediction_quantiles[1,,index], lty = "dashed", col =2)
lines( summary(pred.2)$prediction_quantiles[3,,index], lty = "dashed", col =2)



