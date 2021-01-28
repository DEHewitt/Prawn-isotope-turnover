# Script to conduct a non-linear least squares fit to isotope turnover data
# Compiled by Matthew D. Taylor, NSW Department of Primary Industries
# Last altered: 04/05/2018

library(ggplot2)
library(car)
library(Matching)
library(doBy)
library(nlstools)

rm(list=ls())
while (dev.cur()>1) dev.off()

setwd("E:/Honours/Data")  
infile <- "METMAC_turnover.csv"   #save a .csv file under this name
Data_table = read.csv(infile , header=TRUE)   
summary(Data_table)

Project_filename <- "METMAC_turnover_Output_SUMMARY_d15N.txt"

# Define function
expFun <- function(Time,Sf, tc, S0) {
  Sf+(S0-Sf)*exp(-tc*Time)
}

## Fit model
METMAC_turnover <-nls(d15N~expFun(Time,Sf,tc,S0), data=Data_table, start = list(Sf = 11, tc = 0.03, S0 = 8))
s <- summary(METMAC_turnover)

title <- "##### NLS fit of parameters"
capture.output(title, file = Project_filename, append = FALSE)
capture.output(s, file = Project_filename, append = TRUE)
blank_line <- " "
capture.output(blank_line, file = Project_filename, append = TRUE)
capture.output(blank_line, file = Project_filename, append = TRUE)

# METMAC_turnover_cofdet <- cor(expFun(Data_table$Time,Sf=coef(METMAC_turnover),tc=coef(METMAC_turnover),S0=coef(METMAC_turnover)),Data_table$d15N)^2
# METMAC_turnover_cofdetrnd <-round(METMAC_turnover_cofdet,digits=2)

METMAC_turnover_confint <- nlsBoot(METMAC_turnover)
summary(METMAC_turnover_confint)
METMAC_turnover_est <- cbind(estimates=coef(METMAC_turnover),confint(METMAC_turnover, level = 0.90))

title <- "##### Bootstrapped confidence intervals"
capture.output(title, file = Project_filename, append = TRUE)
capture.output(METMAC_turnover_est, file = Project_filename, append = TRUE)


# Derive confidence intervals for plotting
Time_pred <- seq(0,90,length.out=90)
METMAC_d15N_pred <- expFun(Time_pred,Sf=coef(METMAC_turnover),tc=coef(METMAC_turnover),S0=coef(METMAC_turnover))
ci5_METMAC_d15N <- ci95_METMAC_d15N <- numeric(length(Time_pred))  

# predict from time 0-90 in steps of 1 for each bootstrap of Sf, tc, S0
coefs <- METMAC_turnover_confint$coefboot

# a 1150 x 956 matrix to hold the predictions 
bootfits <- matrix(NA,nrow=length(Time_pred),ncol=nrow(coefs))    #  initially just full of NA's

# then predict a y value for each bootstrap coefficent
for(i in 1:nrow(coefs))
  bootfits[,i] <- expFun(Time_pred,Sf=coefs[i,'Sf'],tc=coefs[i,'tc'],S0=coefs[i,'S0'])

# 5 and 95% quantile for each time
bootci <- apply(bootfits,1,function(x) quantile(x,c(.05,.5,.95)))
ci5_METMAC_d15N <- bootci[1,]
ci95_METMAC_d15N <- bootci[3,]
PredMETMAC_d15N <- bootci[2,]

## PLOT THE RESULTS FROM ALL CURVES
jpeg("METMAC_turnover_d15N.jpg", width = 2048, height = 2048, units = "px", res = 300, pointsize = 12)
# par(mfrow=c(1,1))
par(mar = c(5,5,2,2))

xlim = c(0,90)
ylim = c(0,15)
ylab = expression(delta^15*'N')
xlab = "Time (days)"
plot(Time_pred,PredMETMAC_d15N, xlim = xlim, ylim = ylim, ylab = ylab, xlab = xlab, type = "l", col = "black")
par(new=TRUE)
plot(Time_pred,ci5_METMAC_d15N, axes = F, xlim = xlim, ylim = ylim, ylab = "", xlab = "", type = "l", lty = "dashed", col = "red")
par(new=TRUE)
plot(Time_pred,ci95_METMAC_d15N, axes = F, xlim = xlim, ylim = ylim, ylab = "", xlab = "", type = "l", lty = "dashed", col = "red")
par(new=TRUE)
plot(Data_table$Time,Data_table$d15N, axes = F, xlim = xlim, ylim = ylim, ylab = "", xlab = "", col = "black", cex = 1.5)

dev.off()



