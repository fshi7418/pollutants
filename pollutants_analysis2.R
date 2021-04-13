library(glmnet)
library(corrplot)
pollutants_original <- read.csv("pollutants.csv", header = TRUE)

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


gender_list <- c('female', 'male')
edu_list <- c('before_high_school', 'high_school', 'college', 'college_grad')
eth_list <- c('other', 'mexi_us', 'nonhisp_black', 'nonhisp_white')
smoke_list <- c('no', 'yes')

gender <- gender_list[pollutants$male + 1]
education <- edu_list[pollutants$edu_cat]
race <- eth_list[pollutants$race_cat]
smoke_now <- smoke_list[pollutants$smokenow + 1]

pollutants$male <- factor(gender)
pollutants$edu_cat <- factor(education)
pollutants$race_cat <- factor(race)
pollutants$smokenow <- factor(smoke_now)

pollutants$male <- factor(gender, levels=gender_list)
pollutants$edu_cat <- factor(education, levels=edu_list)
pollutants$race_cat <- factor(race, levels=eth_list)
pollutants$smokenow <- factor(smoke_now, levels=smoke_list)

# basic data manipulation complete


#====
# test the methods
kfold <- 12
pollutants_original$index <- rep(1:kfold, each=N/kfold)

# stepwise selection using AIC
stepAIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants_original$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=pollutants[train.ind, ])
  i_MstepAIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "both", k = 2)
  print(row.names(summary(i_MstepAIC)$coefficients))
  
  i_newdata <- pollutants[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MstepAIC.res <- mspe(i_MstepAIC, i_newdata, i_y)
  
  stepAIC.mspe[i] <- i_MstepAIC.res
}


# stepwise selection using BIC
stepBIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants_original$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=pollutants[train.ind, ])
  i_MstepBIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "both", k = log(length(train.ind)))
  print(row.names(summary(i_MstepBIC)$coefficients))
  
  i_newdata <- pollutants[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MstepBIC.res <- mspe(i_MstepBIC, i_newdata, i_y)
  
  stepBIC.mspe[i] <- i_MstepBIC.res
}


covariate_set1 <- c('ageyrs', 'POP_furan3', 'male', 'ln_lbxcot', 'edu_cat', 'monocyte_pct', 'BMI', 'race_cat',
                   'neutrophils_pct')
covariate_set11 <- covariate_set1
covariate_set1 <- c('length', covariate_set1)
pollutants_set1 <- pollutants[colnames(pollutants) %in% covariate_set1]
model1 <- lm(length ~ ., data=pollutants_set1)

pollutants_set11 <- pollutants[covariate_set11]
vif_list11 <- c()
for (i in 1:length(covariate_set11)) {
  print(covariate_set11[i])
  i_model <- lm(as.formula(paste0(covariate_set11[i], "~ .")), data=pollutants_set11)
  i_vif <- 1 / (1 - summary(i_model)$r.squared)
  vif_list11 <- c(vif_list11, i_vif)
  print(paste0("VIF=", i_vif))
}

covariate_set2 <- c('ageyrs', 'POP_furan3', 'male', 'ln_lbxcot', 'monocyte_pct', 'BMI', 'neutrophils_pct')
covariate_set2 <- c('length', covariate_set2)
pollutants_set2 <- pollutants[colnames(pollutants) %in% covariate_set2]
model2 <- lm(length ~ ., data=pollutants_set2)

# residual analysis
reslin <- resid(model2) # raw residuals
studlin <- reslin/(sigma(model2)*sqrt(1-hatvalues(model2))) # studentized residuals

## plot of residuals vs X
plot(reslin~pollutants_set2$BMI,
     xlab=expression('BMI'),
     ylab="Residuals",
     main="Residuals vs BMI")

# plot of residuals against fitted values
plot(studlin~fitted(model2),
     xlab="fitted values",
     ylab="Residuals",
     main="Residuals vs Fitted Values")

## plot of residuals vs X^2
plot(reslin~I(pollutants_set2$BMI^2),
     xlab=expression("BMI^2"),
     ylab="Residuals",
     main="Residuals vs BMI^2")



# adhoc checks
model_age_v_edu <- lm(ageyrs ~ edu_cat, data=pollutants_set1)
summary(model_age_v_edu)

model_age_v_race <- lm(ageyrs ~ race_cat, data=pollutants_set1)
summary(model_age_v_race)

model_male_v_race <- lm(ageyrs ~ male, data=pollutants_set1)
summary(model_male_v_race)

model_ageyrs_v_everything <- lm(ageyrs ~ edu_cat + race_cat + male, data=pollutants_set1)
summary(model_ageyrs_v_everything)
print(1 / (1 - summary(model_ageyrs_v_everything)$r.squared))
vif(model_ageyrs_v_everything)

summary(lm(edu_cat ~ ageyrs + race_cat + male, data=pollutants_original))
