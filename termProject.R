# SETUP ---------------------------------

library('dplyr')
library('chron')
library('ggplot2')
library("depmixS4")
library('stats')
library(devtools)
require(ggbiplot)

getwd()
setwd("E:/Desktop/Homework/New Classes/CMPT 318/Term Project")

# PCA ---------------------------------

DataDf <- read.table("TermProjectData.txt",header = T,sep = ",")

start <- as.Date("1-01-2007",format="%d-%m-%Y")
end   <- as.Date("31-12-2007",format="%d-%m-%Y")

theDate <- start
i = 1
DataDf$Date <- as.Date(DataDf$Date,format="%d/%m/%Y",tz = "UTC")
dfFinal <- data.frame(Date=as.Date(character()),
                      time=as.chron(character()),
                      GAP=character(), 
                      GRP=character(),
                      GI=character(),
                      stringsAsFactors=FALSE) 
while (theDate < end)
{
  Day <- DataDf[DataDf$Date == theDate,]
  Day_window <- Day %>%
    filter((Time > "18:00:00" & Time < "19:00:01")) %>%
    group_by(Time)
  dfFinal <- rbind(dfFinal,Day_window)
  
  theDate <- theDate + 7 
  i = i+1
}


dfSub <- subset(dfFinal,select = -c(Date,Time))
dfSub<-scale(dfSub)

pca <- prcomp((dfSub),retx = true(), scale. = TRUE, center = TRUE)
ggbiplot(pca)
summary(pca)

ggbiplot(pca)

pca.var <- pca$sdev^2

pca.var.per <- round(pca.var/sum(pca.var)*100, 1)








loading_scores <- pca$rotation[,1]
variable_scores <- abs(loading_scores) 
variable_score_ranked <- sort(variable_scores, decreasing=TRUE)
top_10_variables <- names(variable_score_ranked[1:7])

top_10_variables

pca$rotation[top_10_variables,1]

# MARKOV ---------------------------------
{# MARKOV TRAINING ---------------------------------

dftrain <- subset(dfFinal, select = c(Global_active_power,Global_intensity,Voltage,Date,Time))
numtimes <- rep(60,52)
mod1 <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1),data = dftrain, family=list(gaussian(),gaussian()), nstates = 4, ntimes = numtimes)

##modNew <- depmix(list(dftrain$Global_active_power ~ 1,dftest$Global_intensity ~ 1),family=list(gaussian(),multinomial("identity")),nstates= 4, ntimes = numtimes)
fm1 <- fit(mod1)

summary(fm1)
print(fm1)


mod <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1), data = dftrain,family=list(gaussian(),gaussian()), nstates = 6, ntimes = numtimes)
fm2 <- fit(mod)

summary(fm2)
print(fm2)


mod <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1), data = dftrain,family=list(gaussian(),gaussian()), nstates = 12, ntimes = numtimes)
fm3 <- fit(mod)

summary(fm3)
print(fm3)


mod <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1), data = dftrain,family=list(gaussian(),gaussian()), nstates = 14, ntimes = numtimes)
fm4 <- fit(mod)

summary(fm4)
print(fm4)


mod <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1), data = dftrain,family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimes)
fm5 <- fit(mod)

summary(fm5)
print(fm5)

mod <- depmix(list(dftrain$Global_active_power ~ 1,dftrain$Global_intensity ~ 1), data = dftrain,family=list(gaussian(),gaussian()), nstates = 18, ntimes = numtimes)
fm6 <- fit(mod)

summary(fm6)
print(fm6)

plot(x = c(BIC(fm1),BIC(fm2),BIC(fm3),BIC(fm4),BIC(fm5),BIC(fm6)), y = c(logLik(fm1),logLik(fm2),logLik(fm3),logLik(fm4),logLik(fm5),logLik(fm6)),col = 1:6,pch = 19,xlab =  "BIC",ylab = "LogLikelihood",type = "b" )+
  legend("bottomright", legend = c(nstates(fm1),nstates(fm2),nstates(fm3),nstates(fm4),nstates(fm5),nstates(fm6)), col = 1:6, pch = 19)
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"))
}

{# MARKOV TESTING ---------------------------------

# I don't think the following 2 lines are needed
start <- as.Date("1-01-2008",format="%d-%m-%Y")
end   <- as.Date("24-12-2009",format="%d-%m-%Y")

theDate <- start
i = 1

dftest <- data.frame(Date=as.Date(character()),
                      time=as.chron(character()),
                      GAP=character(), 
                      GRP=character(),
                      GI=character(),
                      stringsAsFactors=FALSE) 
while (theDate < end)
{
  Day <- DataDf[DataDf$Date == theDate,]
  Day_window <- Day %>%
    filter((Time > "18:00:00" & Time < "19:00:01")) %>%
    group_by(Time)
  dftest <- rbind(dftest,Day_window)
  
  theDate <- theDate + 7 # use every Monday
  i = i+1
}

dftest <- subset(dftest, select = c(Global_active_power,Global_intensity,Date,Time))

numtimes2 <- rep(60,100)

testmod3 <- depmix(list(dftest$Global_active_power ~ 1,dftest$Global_intensity ~ 1),family=list(gaussian(),gaussian()), nstates = 12, ntimes = numtimes2)
# get parameters from estimated model



testmod3 <- setpars(testmod3,getpars(fm3))
# check the state sequence and probabilities
viterbi(testmod3)

fb3 <- forwardbackward(testmod3)
all.equal(-sum(log(fb3$sca)),fb3$logLike)

testmod4 <- depmix(list(dftest$Global_active_power ~ 1,dftest$Global_intensity ~ 1),data = dftest,family=list(gaussian(),gaussian()), nstates = 14, ntimes = numtimes2)
# get parameters from estimated model



testmod4 <- setpars(testmod4,getpars(fm4))
# check the state sequence and probabilities
viterbi(testmod4)
fb4 <- forwardbackward(testmod4)
all.equal(-sum(log(fb4$sca)),fb4$logLike)

testmod5 <- depmix(list(dftest$Global_active_power ~ 1,dftest$Global_intensity ~ 1),data = dftest,family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimes2)
# get parameters from estimated model



testmod5 <- setpars(testmod5,getpars(fm5))
# check the state sequence and probabilities
viterbi(testmod5)

fb5 <- forwardbackward(testmod5)
all.equal(-sum(log(fb5$sca)),fb5$logLike)

testmod6 <- depmix(list(dftest$Global_active_power ~ 1,dftest$Global_intensity ~ 1),data = dftest,family=list(gaussian(),gaussian()), nstates = 18, ntimes = numtimes2)
# get parameters from estimated model



testmod6 <- setpars(testmod6,getpars(fm6))
# check the state sequence and probabilities
viterbi(testmod6)
fb6 <- forwardbackward(testmod6)
all.equal(-sum(log(fb6$sca)),fb6$logLike)

fb3$logLike <- fb3$logLike*(52/100)
fb4$logLike <- fb4$logLike*(52/100)
fb5$logLike <- fb5$logLike*(52/100)
fb6$logLike <- fb6$logLike*(52/100)

plot(x = c(nstates(fm3),nstates(fm4),nstates(fm5),nstates(fm6),nstates(fm3),nstates(fm4),nstates(fm5),nstates(fm6)), y = c(logLik(fm3),logLik(fm4),logLik(fm5),logLik(fm6),fb3$logLike,fb4$logLike,fb5$logLike,fb6$logLike),col = c(1,1,1,1,2,2,2,2),pch = 19,xlab =  "Number of states",ylab = "LogLikelihood")+
  legend("bottom", legend = c("training", "testing"),col = c(1,2), pch = 19)
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"))
}


# ANOMALIES ---------------------------------

{#* ANOMALY SETUP ---------------------------------

AnomData1 <- read.table("DataWithAnomalies1.txt",header = T,sep = ",")
AnomData2 <- read.table("DataWithAnomalies2.txt",header = T,sep = ",")
AnomData3 <- read.table("DataWithAnomalies3.txt",header = T,sep = ",")

AnomData1$Date <- as.Date(AnomData1$Date,format="%d/%m/%Y",tz = "UTC")
AnomData2$Date <- as.Date(AnomData2$Date,format="%d/%m/%Y",tz = "UTC")
AnomData3$Date <- as.Date(AnomData3$Date,format="%d/%m/%Y",tz = "UTC")

start <- as.Date("04-01-2010",format="%d-%m-%Y")
end   <- as.Date("27-11-2010",format="%d-%m-%Y")

numtimes3 <- rep(60,47)
numtimesDaily <- 60

}

{# ANOMALOUS DATASET 1 ---------------------------------

theDate <- start
i = 1

anomTest <- data.frame(Date=as.Date(character()),
                     time=as.chron(character()),
                     GAP=character(), 
                     GRP=character(),
                     GI=character(),
                     stringsAsFactors=FALSE) 
while (theDate < end)
{
  Day <- AnomData1[AnomData1$Date == theDate,]
  Day_window <- Day %>%
    filter((Time > "18:00:00" & Time < "19:00:01")) %>%
    group_by(Time)
  anomTest <- rbind(anomTest,Day_window)
  
  theDate <- theDate + 7 # use every Monday
  i = i+1
}

anomTest <- subset(anomTest, select = c(Global_active_power,Global_intensity,Date,Time))

anomTestMod1 <- depmix(list(anomTest$Global_active_power ~ 1,anomTest$Global_intensity ~ 1),family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimes3)
# get parameters from estimated model

anomTestMod1 <- setpars(anomTestMod1,getpars(fm5))  # We use fm5 because it has 16 states
# check the state sequence and probabilities
viterbi(anomTestMod1)

fb5Anom <- forwardbackward(anomTestMod1)
all.equal(-sum(log(fb5Anom$sca)),fb5Anom$logLike)
}
{# ANOMALOUS DATASET 2 ---------------------------------
  
  theDate <- start
  i = 1
  
  anomTest2 <- data.frame(Date=as.Date(character()),
                          time=as.chron(character()),
                          GAP=character(), 
                          GRP=character(),
                          GI=character(),
                          stringsAsFactors=FALSE) 
  
  while (theDate < end)
  {
    Day <- AnomData2[AnomData2$Date == theDate,]
    Day_window <- Day %>%
      filter((Time > "18:00:00" & Time < "19:00:01")) %>%
      group_by(Time)
    anomTest2 <- rbind(anomTest2,Day_window)
    
    theDate <- theDate + 7 # use every Monday
    i = i+1
  }
  
  anomTest2 <- subset(anomTest2, select = c(Global_active_power,Global_intensity,Date,Time))
  
  anomTestMod2 <- depmix(list(anomTest2$Global_active_power ~ 1,anomTest2$Global_intensity ~ 1),family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimes3)
  # get parameters from estimated model
  
  anomTestMod2 <- setpars(anomTestMod2,getpars(fm5))  # We use fm5 because it has 16 states
  # check the state sequence and probabilities
  viterbi(anomTestMod2)
  
  fb5Anom2 <- forwardbackward(anomTestMod2)
  all.equal(-sum(log(fb5Anom2$sca)),fb5Anom2$logLike)
}
{# ANOMALOUS DATASET 3 ---------------------------------

theDate <- start
i = 1

anomTest3 <- data.frame(Date=as.Date(character()),
                       time=as.chron(character()),
                       GAP=character(), 
                       GRP=character(),
                       GI=character(),
                       stringsAsFactors=FALSE) 

while (theDate < end)
{
  Day <- AnomData3[AnomData3$Date == theDate,]
  Day_window <- Day %>%
    filter((Time > "18:00:00" & Time < "19:00:01")) %>%
    group_by(Time)
  anomTest3 <- rbind(anomTest3,Day_window)
  
  theDate <- theDate + 7 # use every Monday
  i = i+1
}

anomTest3 <- subset(anomTest3, select = c(Global_active_power,Global_intensity,Date,Time))

anomTestMod3 <- depmix(list(anomTest3$Global_active_power ~ 1,anomTest3$Global_intensity ~ 1),family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimes3)
# get parameters from estimated model

anomTestMod3 <- setpars(anomTestMod3,getpars(fm5))  # We use fm5 because it has 16 states
# check the state sequence and probabilities
viterbi(anomTestMod3)

fb5Anom3 <- forwardbackward(anomTestMod3)
all.equal(-sum(log(fb5Anom3$sca)),fb5Anom3$logLike)
}

{# DAILY LOG LIKELIHOOD ---------------------------------
  DailyTest <- data.frame(Date=as.Date(character()),
                          time=as.chron(character()),
                          GAP=character(),
                          GRP=character(),
                          GI=character(),
                          stringsAsFactors=FALSE)

  theDate <- start
  i = 1
  while (theDate < end)
  {
    Day <- AnomData3[AnomData3$Date == theDate,]
    Day_window <- Day %>%
      filter((Time > "18:00:00" & Time < "19:00:01")) %>%
      group_by(Time)
    DailyTest <- rbind(DailyTest,Day_window)

    theDate <- theDate + 7 # use every Monday
    i = i+1
  }

  DailyTest <- subset(DailyTest, select = c(Global_active_power,Global_intensity,Date,Time))

  i = 1
  theDate <- start
  logDaily <- rep(0.05, 47)

  while (theDate < end)
  {
    Day <- AnomData3[AnomData3$Date == theDate,]
    Day_window <- Day %>%
      filter((Time > "18:00:00" & Time < "19:00:01")) %>%
      group_by(Time)

    theDate <- theDate + 7 # use every Monday

    DailyTestMod <- depmix(list(Day_window$Global_active_power ~ 1,Day_window$Global_intensity ~ 1),family=list(gaussian(),gaussian()), nstates = 16, ntimes = numtimesDaily)
    # get parameters from estimated model

    DailyTestMod <- setpars(DailyTestMod,getpars(fm5))  # We use fm5 because it has 16 states
    # check the state sequence and probabilities
    viterbi(DailyTestMod)

    fb5Daily <- forwardbackward(DailyTestMod)
    all.equal(-sum(log(fb5Daily$sca)),fb5Daily$logLike)

    logDaily[i] <- fb5Daily$logLike*52

    i = i+1
  }

}



{# PLOT ANOMALIES ---------------------------------
  
# Plot Daily for Set 3:
x <- seq(from = 1, to = 47, by = 1)
plot(x,logDaily, xlim = c(0,50), xlab="Day of Test", ylab="Log Likelihood", main="Daily LL For Set 3", pch=19, col="blue")
grid(col = "lightgray", lty = "dotted", lwd = par("lwd"))

# Plot Entire Sets:
plot(x = c(0, 1, 2, 3), y = c(logLik(fm5),fb5Anom$logLike,fb5Anom2$logLike,fb5Anom3$logLike),col = c(1,2,3,4),pch = 19,xlab =  "Dataset",ylab = "Log_Likelihood")+
  legend("bottomleft", legend = c("Training", "Anomalous1", "Anomalous2", "Anomalous3"),col = c(1,2,3,4), pch = 19)
  grid(col = "lightgray", lty = "dotted", lwd = par("lwd"))
  


  
}

