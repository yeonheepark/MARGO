## Define Simulation Function for parallel processing - currently just take argument for informative or vague prior
simulate <- function(prior_diag,rndmethod,betavec,gammavec,train_control){
  ## Sample outcome based on covariates from Bernoulli distribution
  x1all <- rbinom(N_total, 1, px1)
  x2all <- rbinom(N_total, 1, 0.5)
  
  ## Split sampled data into different cohort sizes
  x1 <- x1all[1:n1]
  x1.2 <- x1all[(n1 + 1):(n1 + n2)]
  x1.3 <- x1all[(n1 + n2 + 1):N_total]
  x2 <- x2all[1:n1]
  x2.2 <- x2all[(n1 + 1):(n1 + n2)]
  x2.3 <- x2all[(n1 + n2 + 1):N_total]
  
  ## Define bound index for interim splits
  bound_index = 1
  
  ## Sample treatment indicator from Bernoulli (1 for A, 0 for B)
  G1 <- rbinom(n1, 1, 0.5)
  
  ## Define vector of 1's for intercept coefficient for probit model (R)
  onevec <- rep(1, n1)
  
  ## Define data for probability sampling based on probit model of treatment outcomes
  xmat1 <- cbind(onevec, x1, x2)
  pn <- dim(xmat1)[2]
  p1 <- pnorm(xmat1 %*% betavec + G1*xmat1 %*% gammavec)
  
  ## Sample binary response of success or failure of treatment
  y1 <- rbinom(n1, 1, p1)
  
  ## Define data based on covariate and treatment indicator
  xdata1 <- cbind(onevec, x1, x2, G1, G1*x1, G1*x2)
  
  ## Fit probit model based on covariate and treatment against the success or failure of treatment
  fit <- glm(y1 ~ xdata1-1, family = binomial(link = probit))
  
  ## Define theta and do Gibb's sampling based on model coefficients with informative or vague prior
  mle_theta <- fit$coefficients
  priorb = list(beta = rep(0, length(mle_theta)), 
                P = diag(rep(prior_diag, length(mle_theta))))
  res1 = bayes.probit(y1, xdata1, nsample, priorb)
  
  #### Gibb's sampling ####
  resbetag <- res1$beta[-(1:nburn),]
  resbeta <- res1$beta[-(1:nburn),1:pn]
  resbetahat <- apply(resbetag, 2, mean)
  
  ## Define difference of patients for Treatment A v B
  ndiff1 <- length(which(G1 == 1)) - length(which(G1 == 0)) 
  
  ## Define difference of failure rate
  nf1 <- length(y1) - sum(y1) 
  
  ## Create data frame for treatment A v B and outcome
  data_total = data <- data.frame()
  data <- data.frame(treatment = G1, outcome = y1)
  data_total <- rbind(data_total, data)
  data_total$treatment <- as.factor(data_total$treatment)
  data_total$outcome <- as.factor(data_total$outcome)
  
  ## Define bound index within data frame
  data_total <- data_total %>% mutate(time = factor(rep(1:bound_index,group[1:bound_index])))
  
  ## Define proportions of patients designated into control and treatment
  ctrl_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 0])))
  trt_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 1])))
  
  ## Clean data for ML modeling
  ml.data <- data_total %>%
    mutate(
      outcome = factor(outcome,levels = c(0,1), labels = c("Success","Fail")),
      x1 = c(x1),
      x2 = c(x2)
    )
  ml.data <- select(ml.data,-3)
  
  ## Define test for superior or futility for interim 1 based on chi square test if treatment proportion is greater than control proportion
  ### and set p value to 1 if not
  if (all(data_total$time == 1) | N_total/block_number < 2) {
    if (((ctrl_prop - trt_prop >= 0) & alternative == "less") | ((trt_prop - ctrl_prop >= 0) & alternative == "greater")) {
      p.val1 <- chisq.test(data_total$treatment, data_total$outcome,correct = correct)$p.value/2
      test1 <- sqrt(as.numeric(chisq.test(data_total$treatment,data_total$outcome, correct = correct)$statistic))
    }
    else {
      p.val1 <- 1
      test1 <- 0
    }
  }else {
    p.val1 <- mantelhaen.test(table(data_total), alternative = alternative,correct = correct)$p.val
    test1 <- sqrt(as.numeric(mantelhaen.test(table(data_total), 
                                             alternative = alternative, correct = correct)$statistic))
  }
  
  ## Perform early stopping test for superiority of futility based on O'brien-fleming alpha spending function
  if (test1 > supbound[bound_index]) {
    ind <- bound_index
    ndiff <- ndiff1
    nf <- nf1
    ind.power <- 1
    return(c(ind,ind.power,ndiff,nf))
    next
  }else if (test1 < futbound[bound_index]) {
    ind <- bound_index 
    ndiff <- ndiff1
    nf <- nf1
    ind.power <- 0
    return(c(ind,ind.power,ndiff,nf))
    next
  }
  
  ## Define next cohort for cohort 2 with bound index 2 if passed early stopping
  else{
    bound_index <- 2
    onevec.2 <- rep(1, n2)
    newxdata2 <- cbind(onevec.2, x1.2, x2.2, onevec.2, x1.2, x2.2)
    newxdatas2 <- newxdata2[,1:pn]
    
    ## Define Cohort Covariate data for ML Modeling
    c2data <- data.frame(x1 = x1.2, x2 = x2.2)
    
    ## Define G2 for treatment indicator for A or B
    G2 <- c()
    
    ## Randomize patients into treatment categories based on probit model using initial cohort data
    ### and the different CARA methods
    if (rndmethod == "SVM") {
      svm_grid <- expand.grid(C = c(2^(-5:5)), sigma = c(2^(-15:-7)))
      alloc_model <- train(outcome ~., data = ml.data, 
                           method = "svmRadial",
                           metric = "ROC",
                           trControl = train_control,
                           tuneGrid = svm_grid)
    }
    if (rndmethod == "KNN") {
      alloc_model <- train(outcome ~., data = ml.data, 
                           method = "knn",
                           metric = "ROC",
                           trControl = train_control,
                           tuneLength = 10)
    }
    if (rndmethod == "RF") {
      rf_grid = expand.grid(mtry = c(1:3))
      alloc_model <- train(outcome ~., data = ml.data, 
                           method = "rf",
                           metric = "ROC",
                           trControl = train_control,
                           tuneGrid = rf_grid)
    }
    if (rndmethod == "ADA") {
      ada_grid <- expand.grid(iter = seq(200,350,50),
                              maxdepth = 1,
                              nu = 0.01)
      alloc_model <- train(outcome ~., data = ml.data, 
                           method = "ada",
                           metric = "ROC",
                           trControl = train_control,
                           tuneGrid = ada_grid)
    }
    if (rndmethod == "NN") {
      nn_grid <- expand.grid(size = 1:10)
      
      alloc_model <- train(outcome ~., data = ml.data, 
                           method = "mlp",
                           metric = "ROC",
                           trControl = train_control,
                           tuneGrid = nn_grid)
    }
    for (l in 1:n2) {
      pred <- c2data[l,]
      pred$treatment <- 0
      pred$treatment <- as.factor(pred$treatment)
      g0 <- as.numeric(predict(alloc_model,newdata = pred))
      
      pred$treatment <- 1
      pred$treatment <- as.factor(pred$treatment)
      g1 <- as.numeric(predict(alloc_model,newdata = pred))
      
      if(g0 > g1){
        G2[l] <- 0
      }else if(g0 < g1){
        G2[l] <- 1
      }else{
        G2[l] <- rbinom(1, 1, 0.5)
      }
    }    
    ## Define data for probability sampling based on probit model of treatment outcomes
    xmat2 <- cbind(onevec.2, x1.2, x2.2)
    p1.2 <- pnorm(xmat2 %*% betavec + G2*xmat2 %*% gammavec)
    
    ## Sample binary response of success or failure of treatment
    y2 <- rbinom(n2, 1, p1.2)
    
    ## Define data based on covariate and treatment indicator
    xdata22 <- cbind(onevec.2, x1.2, x2.2, G2, G2*x1.2, G2*x2.2) 
    
    ## Combine data for treatment and outcome between the 1st and 2nd cohorts
    xdata2 <- rbind(xdata1, xdata22)
    y2a <- c(y1, y2)
    G2a <- c(G1, G2)
    
    ## Fit probit model based on covariates and treatment against the success or failure of treatment
    fit2 <- glm(y2a ~ xdata2-1, family = binomial(link = probit))
    
    ## Define theta and do Gibb's sampling based on model coefficients with informative or vague prior
    mle_theta2 <- fit2$coefficients
    priorb = list(beta = rep(0, length(mle_theta2)), 
                  P = diag(rep(prior_diag, length(mle_theta2))))
    res2 = bayes.probit(y2a, xdata2, nsample, priorb)
    
    #### Gibb's sampling ####
    resbetag2 <- res2$beta[-(1:nburn),]
    resbeta2 <- res2$beta[-(1:nburn),1:pn]
    resbetahat2 <- apply(resbetag2, 2, mean)
    
    ## Define difference of patients for Treatment A v B
    ndiff2 <- length(which(G2 == 1)) - length(which(G2 == 0)) 
    
    ## Define difference of failure rate
    nf2 <- length(y2) - sum(y2)
    data <- data.frame(treatment = G2, outcome = y2)
    
    ## Create data frame for treatment A v B and outcome
    data_total <- rbind(data_total[,1:2], data)
    data_total$treatment <- as.factor(data_total$treatment)
    data_total$outcome <- as.factor(data_total$outcome)
    
    ## Define bound index within data frame
    data_total <- data_total %>% mutate(time = factor(rep(1:bound_index,group[1:bound_index])))
    
    ## Define proportions of patients designated into control and treatment
    ctrl_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 0])))
    trt_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 1])))
    
    ## Clean data for ML modeling
    ml.data <- data_total %>%
      mutate(
        outcome = factor(outcome,levels = c(0,1), labels = c("Success","Fail")),
        x1 = c(x1,x1.2),
        x2 = c(x2,x2.2)
      )
    ml.data <- select(ml.data,-3)
    
    ## Reformat data for Propensity Score Overlap Weighting
    ps.data <- data_total %>%
      mutate(
        treatment = as.numeric(treatment),
        outcome = as.numeric(outcome) - 1,
        x1 = c(x1,x1.2),
        x2 = c(x2,x2.2)
      )
    
    ## Remove time variable from propensity score data
    ps.data <- select(ps.data,-3)
    
    ## Define formula for treatment based on covariates
    ps.form <- treatment ~ x1 + x2
    
    ## Define test for superior or futility for interim 2 based on propensity score weighting method:
    if (all(data_total$time == 1) | N_total/block_number < 2) {
      if (((ctrl_prop - trt_prop >= 0) & alternative == "less") | ((trt_prop - ctrl_prop >= 0) & alternative == "greater")) {
        p.val2 <- chisq.test(data_total$treatment, data_total$outcome, 
                             correct = correct)$p.value/2
        test2 <- sqrt(as.numeric(chisq.test(data_total$treatment, 
                                            data_total$outcome, correct = correct)$statistic))
      }
      else {
        p.val2 <- 1
        test2 <- 0
      }
    }else {
      ## Define propensity score overlap weight R object based on ps.data
      pswgt <- PSweight(ps.formula = ps.form, weight = "overlap",yname = "outcome",data = ps.data)
      
      ## Summarize propensity score weighting
      ato_sum <- summary(pswgt, CI = FALSE)
      
      ## Define p value and test statistic for early stopping based on propensity score data
      p.val2 <- ato_sum$estimates[6]/2
      test2 <- ato_sum$estimates[3]
    }
    
    ## Perform early stopping test for superiority of futility based on O'brien-fleming alpha spending function
    if (test2 > supbound[bound_index]) {
      ind <- bound_index
      ndiff <- ndiff1 + ndiff2
      nf <- nf1 + nf2
      ind.power <- 1
      return(c(ind,ind.power,ndiff,nf))
      next
    }else if (test2 < futbound[bound_index]) {
      ind <- bound_index 
      ndiff <- ndiff1 + ndiff2
      nf <- nf1 + nf2
      ind.power <- 0
      return(c(ind,ind.power,ndiff,nf))
      next
    }
    ## Define next cohort for cohort 3 with bound index 3 if passed early stopping
    else{
      bound_index <- 3
      onevec.3 <- rep(1, n3)
      newxdata3 <- cbind(onevec.3, x1.3, x2.3, onevec.3, x1.3, x2.3)
      newxdatas3 <- newxdata3[,1:pn]
      
      ## Define Cohort Covariate data for ML Modeling
      c3data <- data.frame(x1 = x1.3, x2 = x2.3)
      
      ## Define G3 for treatment indicator for A or B
      G3 <- c()
      
      ## Randomize patients into treatment categories based on probit model using initial cohort data
      ### and the different CARA methods
      if (rndmethod == "SVM") {
        svm_grid <- expand.grid(C = c(2^(-5:5)), sigma = c(2^(-15:-7)))
        alloc_model <- train(outcome ~., data = ml.data, 
                             method = "svmRadial",
                             metric = "ROC",
                             trControl = train_control,
                             tuneGrid = svm_grid)
      }
      if (rndmethod == "KNN") {
        alloc_model <- train(outcome ~., data = ml.data, 
                             method = "knn",
                             metric = "ROC",
                             trControl = train_control,
                             tuneLength = 10)
      }
      if (rndmethod == "RF") {
        rf_grid = expand.grid(mtry = c(1:3))
        alloc_model <- train(outcome ~., data = ml.data, 
                             method = "rf",
                             metric = "ROC",
                             trControl = train_control,
                             tuneGrid = rf_grid)
      }
      if (rndmethod == "ADA") {
        ada_grid <- expand.grid(iter = seq(200,350,50),
                                maxdepth = 1,
                                nu = 0.01)
        alloc_model <- train(outcome ~., data = ml.data, 
                             method = "ada",
                             metric = "ROC",
                             trControl = train_control,
                             tuneGrid = ada_grid)
      }
      if (rndmethod == "NN") {
        nn_grid <- expand.grid(size = 1:10)
        
        alloc_model <- train(outcome ~., data = ml.data, 
                             method = "mlp",
                             metric = "ROC",
                             trControl = train_control,
                             tuneGrid = nn_grid)
      }
      for (l in 1:n3) {
        pred <- c3data[l,]
        pred$treatment <- 0
        pred$treatment <- as.factor(pred$treatment)
        g0 <- as.numeric(predict(alloc_model,newdata = pred))
        
        pred$treatment <- 1
        pred$treatment <- as.factor(pred$treatment)
        g1 <- as.numeric(predict(alloc_model,newdata = pred))
        
        if(g0 > g1){
          G3[l] <- 0
        }else if(g0 < g1){
          G3[l] <- 1
        }else{
          G3[l] <- rbinom(1, 1, 0.5)
        }
      }    
      ## Define data for probability sampling based on probit model of treatment outcomes
      xmat3 <- cbind(onevec.3, x1.3, x2.3)
      p1.3 <- pnorm(xmat3 %*% betavec + G3*xmat3 %*% gammavec)
      
      ## Sample binary response of success or failure of treatment
      y3 <- rbinom(n3, 1, p1.3)
      
      ## Define data based on covariate and treatment indicator
      xdata33 <- cbind(onevec.3, x1.3, x2.3, G3, G3*x1.3, G3*x2.3) 
      xdata33ab <- cbind(onevec.3, x1.3, x2.3, onevec.3, x1.3, x2.3) 
      
      ## Combine data for treatment and outcome between the 1st and 2nd cohorts
      xdata3 <- rbind(xdata1, xdata22, xdata33)
      y3a <- c(y1, y2, y3)
      G3a <- c(G1, G2, G3)
      
      ## Fit probit model based on covariates and treatment against the success or failure of treatment
      fit3 <- glm(y3a ~ xdata3-1, family = binomial(link = probit))
      
      ## Define theta and do Gibb's sampling based on model coefficients with informative or vague prior
      mle_theta3 <- fit3$coefficients
      priorb = list(beta = rep(0, length(mle_theta3 )), 
                    P = diag(rep(prior_diag, length(mle_theta3 ))))
      res3 = bayes.probit(y3a, xdata3, nsample, priorb)
      
      #### Gibb's sampling ####
      resbetag3 <- res3$beta[-(1:nburn),]
      resbeta3 <- res3$beta[-(1:nburn),1:pn]
      
      ## Define difference of patients for Treatment A v B
      ndiff3 <- length(which(G3 == 1)) - length(which(G3 == 0)) 
      
      ## Define difference of failure rate
      nf3 <- length(y3) - sum(y3)    
      data <- data.frame(treatment = G3, outcome = y3)
      
      ## Create data frame for treatment A v B and outcome
      data_total <- rbind(data_total[,1:2], data)
      data_total$treatment <- as.factor(data_total$treatment)
      data_total$outcome <- as.factor(data_total$outcome)
      
      ## Define bound index within data frame
      data_total <- data_total %>% mutate(time = factor(rep(1:bound_index,group[1:bound_index])))
      
      ## Define proportions of patients designated into control and treatment
      ctrl_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 0])))
      trt_prop <- mean(as.numeric(as.character(data_total$outcome[data_total$treatment == 1])))
      
      ## Reformat data for Propensity Score Overlap Weighting
      ps.data <- data_total %>%
        mutate(
          treatment = as.numeric(treatment),
          outcome = as.numeric(outcome) - 1,
          x1 = c(x1,x1.2,x1.3),
          x2 = c(x2,x2.2,x2.3)
        )
      
      ## Remove time variable from propensity score data
      ps.data <-   select(ps.data,-3)
      
      ## Define formula for treatment based on covariates
      ps.form <- treatment ~ x1 + x2
      
      ## Define test for superior or futility for interim 2 based on propensity score weighting method:
      if (all(data_total$time == 1) | N_total/block_number < 2) {
        if (((ctrl_prop - trt_prop >= 0) & alternative == "less") | ((trt_prop - ctrl_prop >= 0) & alternative == "greater")) {
          p.val3 <- chisq.test(data_total$treatment, data_total$outcome, 
                               correct = correct)$p.value/2
          test3 <- sqrt(as.numeric(chisq.test(data_total$treatment,data_total$outcome, correct = correct)$statistic))
        }
        else {
          p.val3 <- 1
          test3 <- 0
        }
      }else {
        ## Define propensity score overlap weight R object based on ps.data
        pswgt <- PSweight(ps.formula = ps.form, weight = "overlap",yname = "outcome",data = ps.data)
        
        ## Summarize propensity score weighting
        ato_sum <- summary(pswgt, CI = FALSE)
        
        ## Define p value and test statistic for early stopping based on propensity score data
        p.val3 <- ato_sum$estimates[6]/2
        test3 <- ato_sum$estimates[3]
      }
      if (test3 > supbound[bound_index]) {
        ind <- 3
        ndiff <- ndiff1 + ndiff2 + ndiff3
        nf <- nf1 + nf2 + nf3
        ind.power <- 1
      }else{
        ind <- 3
        ndiff <- ndiff1 + ndiff2 + ndiff3
        nf <- nf1 + nf2 + nf3
        ind.power <- 0
      }
    }
  }
  ### Return values of interim indicator, power, allocation difference, number of failures, p value and the test statistic
  return(c(ind,ind.power,ndiff,nf))
}
