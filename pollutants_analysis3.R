library(glmnet)
library(corrplot)
library(car)
pollutants_original <- read.csv("pollutants.csv", header = TRUE)
N <- nrow(pollutants_original)


#==== explorotary data analysis
num_smokenow <- length(which(pollutants_original$smokenow == 1))
num_nosmokenow <- nrow(pollutants_original) - num_smokenow
print(num_smokenow / N)

num_smoked <- length(which(pollutants_original$yrssmoke > 0))
num_nosmoked <- N - num_smoked
print(num_smoked / N)

# drop the first column, which is just indices
pollutants <- pollutants_original[-1]

# drop a white cell pct due to perfect linearity
pollutants <- pollutants[!colnames(pollutants) %in% c('eosinophils_pct')]

# standardize all POP covariates
for (i in 1:11) {
  i_colname <- paste0('POP_PCB', i)
  pollutants[i_colname] <- (pollutants[[i_colname]] - mean(pollutants[[i_colname]])) / sd(pollutants[[i_colname]])
}
for (i in 1:3) {
  i_colname <- paste0('POP_dioxin', i)
  pollutants[i_colname] <- (pollutants[[i_colname]] - mean(pollutants[[i_colname]])) / sd(pollutants[[i_colname]])
}
for (i in 1:4) {
  i_colname <- paste0('POP_furan', i)
  pollutants[i_colname] <- (pollutants[[i_colname]] - mean(pollutants[[i_colname]])) / sd(pollutants[[i_colname]])
}

# deal with categorical covariates
gender_list <- c('female', 'male')
edu_list <- c('before_high_school', 'high_school', 'college', 'college_grad')
eth_list <- c('other', 'mexi_us', 'nonhisp_black', 'nonhisp_white')
smoke_list <- c('no', 'yes')

gender <- gender_list[pollutants$male + 1]
education <- edu_list[pollutants$edu_cat]
race <- eth_list[pollutants$race_cat]
smoke_now <- smoke_list[pollutants$smokenow + 1]

pollutants$male <- factor(gender, levels=gender_list)
pollutants$edu_cat <- factor(education, levels=edu_list)
pollutants$race_cat <- factor(race, levels=eth_list)
pollutants$smokenow <- factor(smoke_now, levels=smoke_list)

# basic data manipulation complete


#==== train test split
set.seed(57)
# 70% - 30% split
sampleTrain <- sample(1:N, round(N*0.7,0), replace = FALSE)

dataTrain <- pollutants[sampleTrain,]
dataTest <- pollutants[-sampleTrain,]


#====
# test the selection algorithms using cross validation (steps 5 and 6 in section 4 of the report)
Ntrain <- nrow(dataTrain)
kfold <- 11  # 11 because Ntrain is divisible by 11, and 11 is close to the convention of 10-fold CV
dataTrain_original <- dataTrain
dataTrain_original$index <- rep(1:kfold, each=Ntrain/kfold)

# forward selection using AIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

# keep track of the MSPE of each fold and storing them
mspe_fAIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ]) # empty model
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ]) # full model with no selection
  i_MfwdAIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "forward", k = 2)
  i_coef_df <- data.frame(summary(i_MfwdAIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1  # indicator of which covariates are selected
  for (j in 2:dim(i_coef_df)[1]) {  # store the p-values of selected covariates in the matrix
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MfwdAIC)$coefficients))

  # compute the MPSE of the i-th fold
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MfwdAIC, newdata=i_newdata)
  mspe_fAIC[i] <- mean(i_res^2)
}
avgMSPE_fAIC <- mean(mspe_fAIC)
write.csv(data.frame(pval_mat), "fwdAIC_pval.csv")
write.csv(data.frame(selection_mat), "fwdAIC_indicator.csv")


# repeat all of the above for the remaining five model selection algorithms
# backward elimination using AIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_bAIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MbckAIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "back", k = 2)
  i_coef_df <- data.frame(summary(i_MbckAIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MbckAIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MbckAIC, newdata=i_newdata)
  mspe_bAIC[i] <- mean(i_res^2)
}
write.csv(data.frame(pval_mat), "bckAIC_pval.csv")
write.csv(data.frame(selection_mat), "bckAIC_indicator.csv")


# stepwise selection using AIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_sAIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MstepAIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "both", k = 2)
  i_coef_df <- data.frame(summary(i_MstepAIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MstepAIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MstepAIC, newdata=i_newdata)
  mspe_sAIC[i] <- mean(i_res^2)
}
write.csv(data.frame(pval_mat), "stepAIC_pval.csv")
write.csv(data.frame(selection_mat), "stepAIC_indicator.csv")


# stepwise selection using BIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_sBIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MstepBIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "both", k = log(length(train.ind)))
  i_coef_df <- data.frame(summary(i_MstepBIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MstepBIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MstepBIC, newdata=i_newdata)
  mspe_sBIC[i] <- mean(i_res^2)
}
write.csv(data.frame(pval_mat), "stepBIC_pval.csv")
write.csv(data.frame(selection_mat), "stepBIC_indicator.csv")


# backward elimination using BIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_bBIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MbckBIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "backward", k = log(length(train.ind)))
  i_coef_df <- data.frame(summary(i_MbckBIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MbckBIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MbckBIC, newdata=i_newdata)
  mspe_bBIC[i] <- mean(i_res^2)
}
write.csv(data.frame(pval_mat), "bckBIC_pval.csv")
write.csv(data.frame(selection_mat), "bckBIC_indicator.csv")


# forward selection using BIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_fBIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MfwdBIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "forward", k = log(length(train.ind)))
  i_coef_df <- data.frame(summary(i_MfwdBIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MfwdBIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MfwdBIC, newdata=i_newdata)
  mspe_fBIC[i] <- mean(i_res^2)
}
write.csv(data.frame(pval_mat), "fwdBIC_pval.csv")
write.csv(data.frame(selection_mat), "fwdBIC_indicator.csv")


#==== steps 5 and 6 complete

# arriving at the final model
# model3: make sense of the effect of cotinine
covariate_set3 <- c('ageyrs', 'POP_furan3', 'ln_lbxcot', 'yrssmoke', 'POP_PCB1', 'POP_PCB5', 'edu_cat')
covariate_set3 <- c('length', covariate_set3)
pollutants_set3 <- dataTrain
pollutants_set3$yrssmoke <- as.integer(pollutants_set3$yrssmoke > 0)

# make the yrssmoked covariate binary
smoke_list <- c('never_smoked', 'smoked')
smoke_now <- smoke_list[pollutants_set3$yrssmoke + 1]
pollutants_set3$yrssmoke <- factor(smoke_now, levels=smoke_list)
print(dim(pollutants_set3))

# select only a the covariates that were chosen repeatedly in steps 5 and 6
pollutants_set3 <- pollutants_set3[colnames(pollutants_set3) %in% covariate_set3]

# log transforming the response variable
pollutants_set3$length <- log(pollutants_set3$length)
model3 <- lm(length ~ ageyrs + POP_furan3 + ln_lbxcot + edu_cat + yrssmoke + yrssmoke * ln_lbxcot, 
             data=pollutants_set3)
summary(model3)

# proto_model3 shows that ln_lbxcot on its own has no statistical significance
proto_model3 <- lm(length ~ ageyrs + POP_furan3 + ln_lbxcot + edu_cat, data=pollutants_set3)
summary(proto_model3)


# in-sample cross validation of model3
kfold <- 11
model3_mpse <- rep(NA, kfold)
pollutants_set31 <- pollutants_set3
pollutants_set31$index <- rep(1:kfold, each=nrow(pollutants_set3)/kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants_set31$index != i)
  
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_model3 <- update(model3, subset=train.ind)
  
  model3.res <- pollutants_set3$length[-train.ind] - 
    predict(i_model3, newdata = pollutants_set3[-train.ind, ])

  model3_mpse[i] <- mean(model3.res^2) 
}
print(paste0("model3 avg. MPSE:", mean(model3_mpse)))


# residual analysis of model3
reslin <- resid(model3) # raw residuals
studlin <- reslin/(sigma(model3)*sqrt(1-hatvalues(model3))) # studentized residuals

# QQ plot
qqnorm(studlin)
abline(0, 1)
# before log transform of response, the qqplot looked like a gamma distribution


# plot of residuals vs covariates
plotted_col <- "ageyrs"
plot(reslin~pollutants_set3[[plotted_col]],
     xlab=plotted_col,
     ylab="Residuals",
     main=paste0("Residuals vs ", plotted_col))
plotted_col <- "edu_cat"
plot(reslin~pollutants_set3[[plotted_col]],
     xlab=plotted_col,
     ylab="Residuals",
     main=paste0("Residuals vs ", plotted_col))
plotted_col <- "POP_furan3"
plot(reslin~pollutants_set3[[plotted_col]],
     xlab=plotted_col,
     ylab="Residuals",
     main=paste0("Residuals vs ", plotted_col))
plotted_col <- "ln_lbxcot"
plot(reslin~pollutants_set3[[plotted_col]],
     xlab=plotted_col,
     ylab="Residuals",
     main=paste0("Residuals vs ", plotted_col))

# plot of residuals against fitted values
plot(studlin~fitted(model3),
     xlab="fitted values",
     ylab="Residuals",
     main="Residuals vs Fitted Values")


# outliers analysis
# number of covariates
p <- length(model3$coefficients) - 1
pred <- predict(model3, newdata=pollutants_set3)

# check leverage
h <- hatvalues(model2)
which(h > 2* (p + 1) / N)

# DFFITS
dffits_model3 <- dffits(model3) 
## plot DFFITS
plot(dffits_model3,ylab="DFFITS") 
abline(h=2 * sqrt((p + 1) / N), lty=2)  ## add thresholds
abline(h=-2 * sqrt((p + 1) / N),lty=2)
## highlight influential points
dff_ind <- which(abs(dffits_model3) > 2 * sqrt((p + 1) / N))
points(dffits_model3[dff_ind]~dff_ind, col="red", pch=19) ## add red points
# we exclude text because it makes the chart too cluttered
# text(y=dffits_model2[dff_ind], x=dff_ind, labels=dff_ind, pos=2) ## label high influence points


# Cook's Distance

D <- cooks.distance(model3) # Cook's distance
## influential points
inf_ind <- which(pf(D, p + 1, N - p - 1, lower.tail=TRUE) > 0.5)
## plot cook's Distance
plot(D,ylab="Cook's Distance")
points(D[inf_ind]~inf_ind, col="red", pch=19) ## add red points
text(y=D[inf_ind], x=inf_ind, labels=inf_ind, pos=4) ## label high influence points


# DFBETAS
DFBETAS <- dfbetas(model3) 
dim(DFBETAS)
## beta1 (furan3)
plot(DFBETAS[, 2], type="h", xlab="Obs. Number", ylab=expression(paste("DFBETAS: ",beta[1])))
num_outliers <- length(which(abs(DFBETAS[, 2]) > 2 / sqrt(N)))
show_points <- order(-abs(DFBETAS[, 2]))[1:num_outliers] 
points(x=show_points,y=DFBETAS[show_points,2], pch=19, col="red")
text(x=show_points, y=DFBETAS[show_points,2], labels=show_points,pos=2)

## beta2 (gender)
plot(DFBETAS[,3], type="h",xlab="Obs. Number", ylab=expression(paste("DFBETAS: ",beta[2])))
num_outliers <- length(which(abs(DFBETAS[, 3]) > 2 / sqrt(N)))
show_points <- order(-abs(DFBETAS[, 3]))[1:num_outliers] 
points(x=show_points,y=DFBETAS[show_points,3], pch=19, col="red")
text(x=show_points, y=DFBETAS[show_points, 3], labels=show_points, pos=4)

## beta3 (age)
plot(DFBETAS[, 4], type="h",xlab="Obs. Number", ylab=expression(paste("DFBETAS: ", beta[3])))
num_outliers <- length(which(abs(DFBETAS[, 4]) > 2 / sqrt(N)))
show_points <- order(-abs(DFBETAS[, 4]))[1:num_outliers] 
points(x=show_points,y=DFBETAS[show_points, 4], pch=19, col="red")
text(x=show_points, y=DFBETAS[show_points, 4], labels=show_points, pos=4)

## rule of thumb
which(abs(DFBETAS[, 2]) > 2 / sqrt(N))  # furan3
which(abs(DFBETAS[, 3]) > 2 / sqrt(N))  # gender
which(abs(DFBETAS[, 4]) > 2 / sqrt(N))  # ageyrs


#==== we select model3 to apply to the testing set
pollutants_set3_test <- dataTest
pollutants_set3_test$yrssmoke <- as.integer(pollutants_set3_test$yrssmoke > 0)
pollutants_set3_test$length <- log(pollutants_set3_test$length)

smoke_list <- c('never_smoked', 'smoked')
smoke_now <- smoke_list[pollutants_set3_test$yrssmoke + 1]
pollutants_set3_test$yrssmoke <- factor(smoke_now, levels=smoke_list)
print(dim(pollutants_set3_test))

testing_origin <- pollutants_set3_test$length
testing_pred <- predict(model3, newdata=pollutants_set3_test)
model3_test_res <- testing_origin - testing_pred
print(paste0("testing MPSE for model3 is ", mean(model3_test_res^2)))


# baseline models compared to model3
baseline_data <- dataTrain
baseline_data$yrssmoke <- as.integer(baseline_data$yrssmoke > 0)
smoke_list <- c('never_smoked', 'smoked')
smoke_now <- smoke_list[baseline_data$yrssmoke + 1]
baseline_data$yrssmoke <- factor(smoke_now, levels=smoke_list)

baseline_data$length <- log(baseline_data$length)
m0 <- lm(baseline_data$length ~ 1)
mfull <- lm(length ~ ., data=baseline_data)
summary(m0)
summary(mfull)

m0_testing_pred <- predict(m0, newdata=pollutants_set3_test)
m0_test_res <- testing_origin - m0_testing_pred
mfull_testing_pred <- predict(mfull, newdata=pollutants_set3_test)
mfull_test_res <- testing_origin - mfull_testing_pred

print(paste0("testing MPSE for m0 is ", mean(m0_test_res^2)))
print(paste0("testing MPSE for mfull is ", mean(mfull_test_res^2)))


# PCB1, 5 have high correlations with furan3

summary(lm(POP_PCB1 ~ POP_PCB5 + POP_furan3 + edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_PCB1 + POP_furan3 + edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_furan3 ~ edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_PCB1, data=pollutants_original))
summary(lm(POP_PCB1 ~ POP_furan3, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_furan3, data=pollutants_original))





#==== LASSO
datalassoTrain <- dataTrain
datalassoTrain$length <- log(datalassoTrain$length)
datalassoTest <- dataTest
datalassoTest$length <- log(datalassoTest$length)
Mfull <- lm(length ~ ., data=datalassoTrain)
Mfull_test <- lm(length ~ ., data=datalassoTest)
X_train <- model.matrix(Mfull)
y_train <- datalassoTrain[['length']]
X_test <- model.matrix(Mfull_test)
y_test <- datalassoTest[['length']]
M_lasso <- glmnet(x=X_train,y=y_train,alpha = 1)

## plot paths
plot(M_lasso,xvar = "lambda",label=TRUE)

## fit with crossval
cvfit_lasso <-  cv.glmnet(x=X_train,y=y_train,alpha = 1)

## plot MSPEs by lambda
plot(cvfit_lasso)

## estimated betas for minimum lambda 
coef(cvfit_lasso, s = "lambda.min")## alternatively could use "lambda.1se"

## predictions
pred_lasso <- predict(cvfit_lasso,newx=X_test,  s="lambda.min")

## MSPE in test set
MSPE_lasso <- mean((pred_lasso-y_test)^2)
print(MSPE_lasso)
