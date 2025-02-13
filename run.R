### R Code to Run simulations in parallel for MLAL 
## Call libraries
library(ldbounds)
library(tidyverse)
library(LearnBayes)
library(PSweight)
library(caret)
library(doRNG)
library(randomForest)
library(RSNNS)


## Set the seed
set.seed(840130)

## Set cores for parallel processing
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

## Set clinical trial cohort sizes and number of blocks for simulation
N_total <- 210
n1 = n2 = n3 = 70
group <- c(70, 70, 70)
block_number <- 3

## Define group sequential bounds using Lan-DeMets alpha spending function
### Additionally, define superiority and futility bounds in the trial
time <- c(1/3, 2/3, 1)
bounds <- ldBounds(time, iuse = c(1, 1), alpha = c(0.025, 0.025))
supbound <- bounds$upper.bounds 
futbound <- bounds$lower.bounds

## Set response probability of biomarker 1
px1 <- 0.5

## Number of simulations
n_sim <- 1000

## Define number of sampling and burn for Gibb's sampling and testing alternative
nsample = 10000
nburn = 5000
alternative = "greater"
correct = FALSE

### Set training parameters for the models
train_param <- trainControl(method = "cv",
                            number = 5,
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary)

## Set start time for parallel processing run
start_time <- Sys.time()

## Run simulations for different CARAML with OW methods ADA, SVM, KNN, NN and RF 
### and append to final dataframe 

Scenario = "Scenario7"

result_p <- rep(NA,5)
result_diff <- rep(NA,5)
result_nf <- rep(NA,5)
result_estp <- rep(NA,5)

result <- foreach(i = 1:n_sim, .combine = 'cbind',.packages = c('LearnBayes','dplyr','caret','PSweight')) %dorng% simulate(
  prior_diag = 2,
  rndmethod = "ADA",
  betavec = c(0,0,0),
  gammavec = c(0.5,0,0),
  train_control = train_param)

attr(result, 'rng') <- NULL
earlyst1 <- sum(result[1,]==1) + sum(result[1,] == 2)
earlyst2 <- sum(result[1,]==1) + sum(result[1,] == 2) + sum(result[1,] ==3)
earlyst <- earlyst1/earlyst2
result_p[1] <- sum(result[2,],na.rm=T)/earlyst2
result_diff[1] <- mean(result[3,])
result_nf[1] <- mean(result[4,])
result_estp[1] <- earlyst

result <- foreach(i = 1:n_sim, .combine = 'cbind',.packages = c('LearnBayes','dplyr','caret','PSweight')) %dorng% simulate(
  prior_diag = 2,
  rndmethod = "SVM",
  betavec = c(0,0,0),
  gammavec = c(0.5,0,0),
  train_control = train_param)

attr(result, 'rng') <- NULL
earlyst1 <- sum(result[1,]==1) + sum(result[1,] == 2)
earlyst2 <- sum(result[1,]==1) + sum(result[1,] == 2) + sum(result[1,] ==3)
earlyst <- earlyst1/earlyst2
result_p[2] <- sum(result[2,],na.rm=T)/earlyst2
result_diff[2] <- mean(result[3,])
result_nf[2] <- mean(result[4,])
result_estp[2] <- earlyst

result <- foreach(i = 1:n_sim, .combine = 'cbind',.packages = c('LearnBayes','dplyr','caret','PSweight')) %dorng% simulate(
  prior_diag = 2,
  rndmethod = "KNN",
  betavec = c(0,0,0),
  gammavec = c(0.5,0,0),
  train_control = train_param)

attr(result, 'rng') <- NULL
earlyst1 <- sum(result[1,]==1) + sum(result[1,] == 2)
earlyst2 <- sum(result[1,]==1) + sum(result[1,] == 2) + sum(result[1,] ==3)
earlyst <- earlyst1/earlyst2
result_p[3] <- sum(result[2,],na.rm=T)/earlyst2
result_diff[3] <- mean(result[3,])
result_nf[3] <- mean(result[4,])
result_estp[3] <- earlyst

result <- foreach(i = 1:n_sim, .combine = 'cbind',.packages = c('LearnBayes','dplyr','caret','PSweight')) %dorng% simulate(
  prior_diag = 2,
  rndmethod = "RF",
  betavec = c(0,0,0),
  gammavec = c(0.5,0,0),
  train_control = train_param)

attr(result, 'rng') <- NULL
earlyst1 <- sum(result[1,]==1) + sum(result[1,] == 2)
earlyst2 <- sum(result[1,]==1) + sum(result[1,] == 2) + sum(result[1,] ==3)
earlyst <- earlyst1/earlyst2
result_p[4] <- sum(result[2,],na.rm=T)/earlyst2
result_diff[4] <- mean(result[3,])
result_nf[4] <- mean(result[4,])
result_estp[4] <- earlyst

result <- foreach(i = 1:n_sim, .combine = 'cbind',.packages = c('LearnBayes','dplyr','caret','PSweight')) %dorng% simulate(
  prior_diag = 2,
  rndmethod = "NN",
  betavec = c(0,0,0),
  gammavec = c(0.5,0,0),
  train_control = train_param)

attr(result, 'rng') <- NULL
earlyst1 <- sum(result[1,]==1) + sum(result[1,] == 2)
earlyst2 <- sum(result[1,]==1) + sum(result[1,] == 2) + sum(result[1,] ==3)
earlyst <- earlyst1/earlyst2
result_p[5] <- sum(result[2,],na.rm=T)/earlyst2
result_diff[5] <- mean(result[3,])
result_nf[5] <- mean(result[4,])
result_estp[5] <- earlyst

output <- data.frame(Scenario = rep(Scenario,5),ML = c("ADA","SVM","KNN","RF","NN"),RejP = result_p,DiffA = result_diff, NF = result_nf,
                     EarlyStop = result_estp)

output

## StopClusters after done with simulations
parallel::stopCluster(cl = my.cluster)
