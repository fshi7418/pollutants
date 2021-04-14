library(glmnet)
library(corrplot)
pollutants_original <- read.csv("pollutants.csv", header = TRUE)
N <- nrow(pollutants_original)


#==== explorotary data analysis
num_smokenow <- length(which(pollutants_original$smokenow == 1))
num_nosmokenow <- nrow(pollutants_original) - num_smokenow
print(num_smokenow / N)

num_smoked <- length(which(pollutants_original$yrssmoke > 0))
num_nosmoked <- N - num_smoked
print(num_smoked / N)

# drop the first column of indices
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

pollutants$length <- log(pollutants$length)

# basic data manipulation complete


#==== train test split
set.seed(57)
sampleTrain <- sample(1:N, round(N*0.7,0), replace = FALSE)

dataTrain <- pollutants[sampleTrain,]
dataTest <- pollutants[-sampleTrain,]


#====
# test the methods
Ntrain <- nrow(dataTrain)
kfold <- 11
dataTrain_original <- dataTrain
dataTrain_original$index <- rep(1:kfold, each=Ntrain/kfold)

# forward selection using AIC and store them in a matrix
selection_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)  # +6 to rownum bc of categorical covariates
pval_mat <- matrix(0, nrow=length(colnames(dataTrain)) + 6, ncol=kfold)
mat_rownames <- c(colnames(dataTrain), c("edu_cathigh_school", "edu_catcollege", "edu_catcollege_grad", 
                                          "race_catmaxi_us", "race_catnonhisp_black", "race_catnonhisp_white"))
row.names(selection_mat) <- mat_rownames
row.names(pval_mat) <- mat_rownames

mspe_fAIC <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(dataTrain_original$index != i)
  
  M0 <- lm(length ~ 1, data=dataTrain[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=dataTrain[train.ind, ])
  i_MfwdAIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "forward", k = 2)
  i_coef_df <- data.frame(summary(i_MfwdAIC)$coefficients)
  selection_mat[mat_rownames %in% row.names(i_coef_df), i] <- 1
  for (j in 2:dim(i_coef_df)[1]) {
    pval_mat[row.names(i_coef_df)[j] == mat_rownames, i] <- i_coef_df[j, "Pr...t.."]
  }
  print(row.names(summary(i_MfwdAIC)$coefficients))
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_res <- i_y - predict(i_MfwdAIC, newdata=i_newdata)
  mspe_fAIC[i] <- mean(i_res^2)
}
avgMSPE_fAIC <- mean(mspe_fAIC)
write.csv(data.frame(pval_mat), "fwdAIC_pval.csv")
write.csv(data.frame(selection_mat), "fwdAIC_indicator.csv")


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




# first attempt to choose a model and tweaking it
# model1
covariate_set1 <- c('ageyrs', 'POP_furan3', 'male', 'ln_lbxcot', 'edu_cat', 'monocyte_pct', 'BMI', 'race_cat',
                   'neutrophils_pct')
covariate_set11 <- covariate_set1
covariate_set1 <- c('length', covariate_set1)
pollutants_set1 <- dataTrain[colnames(dataTrain) %in% covariate_set1]
model1 <- lm(length ~ ., data=pollutants_set1)

pollutants_set11 <- dataTrain[covariate_set11]
vif_list11 <- c()
for (i in 1:length(covariate_set11)) {
  print(covariate_set11[i])
  i_model <- lm(as.formula(paste0(covariate_set11[i], "~ .")), data=pollutants_set11)
  i_vif <- 1 / (1 - summary(i_model)$r.squared)
  vif_list11 <- c(vif_list11, i_vif)
  print(paste0("VIF=", i_vif))
}


# take out education and race because of their high correlation with ageyrs
# covariate_set2 <- c('ageyrs', 'POP_furan3', 'male', 'ln_lbxcot', 'monocyte_pct', 'BMI', 'neutrophils_pct')
covariate_set2 <- c('ageyrs', 'POP_furan3', 'POP_PCB1', 'POP_PCB5')
covariate_set2 <- c('length', covariate_set2)
pollutants_set2 <- dataTrain[colnames(dataTrain) %in% covariate_set2]
# log transforming the response variable
pollutants_set2$length <- log(pollutants_set2$length)
model2 <- lm(length ~ ., data=pollutants_set2)
model2 <- lm(length ~ ., data=pollutants_set2)
summary(model2)

# cross validate model2 and the bare minimum model with ageyrs + furan3
model_min <- lm(length ~ ageyrs + POP_furan3, data=pollutants_set2)
kfold <- 11
model2_mpse <- rep(NA, kfold)
model_min_mpse <- rep(NA, kfold)
pollutants_set21 <- pollutants_set2
pollutants_set21$index <- rep(1:kfold, each=nrow(pollutants_set2)/kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants_set21$index != i)
  
  i_newdata <- dataTrain[-train.ind, ]
  i_y <- dataTrain$length[-train.ind]
  
  i_model_min <- update(model_min, subset=train.ind)
  i_model2 <- update(model2, subset=train.ind)
  
  model_min.res <- pollutants_set2$length[-train.ind] - 
    predict(i_model_min, newdata = pollutants_set2[-train.ind, ])
  model2.res <- pollutants_set2$length[-train.ind] - 
    predict(i_model2, newdata = pollutants_set2[-train.ind, ])
 
  model_min_mpse[i] <- mean(model_min.res^2)
  model2_mpse[i] <- mean(model2.res^2) 
}
print(paste0("minimum model with only age and furan3 avg. MPSE:", mean(model_min_mpse)))
print(paste0("model2 avg. MPSE:", mean(model2_mpse)))


# model3: make sense of the effect of cotinine
covariate_set3 <- c('ageyrs', 'POP_furan3', 'ln_lbxcot', 'yrssmoke', 'POP_PCB1', 'POP_PCB5', 'edu_cat')
covariate_set3 <- c('length', covariate_set3)
pollutants_set3 <- dataTrain
pollutants_set3$yrssmoke <- as.integer(pollutants_set3$yrssmoke > 0)

smoke_list <- c('never_smoked', 'smoked')
smoke_now <- smoke_list[pollutants_set3$yrssmoke + 1]
pollutants_set3$yrssmoke <- factor(smoke_now, levels=smoke_list)
print(dim(pollutants_set3))
###
pollutants_set3 <- pollutants_set3[colnames(pollutants_set3) %in% covariate_set3]

# log transforming the response variable
pollutants_set3$length <- log(pollutants_set3$length)
model3 <- lm(length ~ ageyrs + POP_furan3 + ln_lbxcot + edu_cat + yrssmoke + yrssmoke * ln_lbxcot, 
             data=pollutants_set3)
summary(model3)

proto_model3 <- lm(length ~ ageyrs + POP_furan3 + ln_lbxcot + edu_cat, data=pollutants_set3)
summary(proto_model3)

# cv model3
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
DFBETAS <- dfbetas(model2) 
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





# adhoc checks
model_age_v_edu <- lm(ageyrs ~ edu_cat, data=pollutants_set1)
summary(model_age_v_edu)

model_length_v_edu <- lm(length ~ edu_cat, data=pollutants_set1)
summary(model_length_v_edu)

model_yrssmoke_v_edu <- lm(yrssmoke ~ edu_cat, data=pollutants)
summary(model_yrssmoke_v_edu)

model_age_v_race <- lm(ageyrs ~ race_cat, data=pollutants_set1)
summary(model_age_v_race)

model_male_v_race <- lm(ageyrs ~ male, data=pollutants_set1)
summary(model_male_v_race)

model_ageyrs_v_everything <- lm(ageyrs ~ edu_cat + race_cat + male, data=pollutants_set1)
summary(model_ageyrs_v_everything)
print(1 / (1 - summary(model_ageyrs_v_everything)$r.squared))
vif(model_ageyrs_v_everything)

summary(lm(edu_cat ~ ageyrs + race_cat + male, data=pollutants_original))


summary(lm(POP_PCB1 ~ POP_PCB5 + POP_furan3 + edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_PCB1 + POP_furan3 + edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_furan3 ~ edu_cat + ln_lbxcot + ageyrs, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_PCB1, data=pollutants_original))
summary(lm(POP_PCB1 ~ POP_furan3, data=pollutants_original))
summary(lm(POP_PCB5 ~ POP_furan3, data=pollutants_original))
