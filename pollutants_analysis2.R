library(glmnet)
library(corrplot)
pollutants_original <- read.csv("pollutants.csv", header = TRUE)
pollutants <- pollutants_original[-1]

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

## Factors removed 
pollutants_factorsRemoved <- pollutants[-c(27,28,29,32)]

pollutants_factorsRemovedCor <- cor(pollutants_factorsRemoved)

set.seed(57)

N <- nrow(pollutants)
sampleTrain <- sample(1:N, round(N*0.8,0), replace = FALSE)

dataTrain <- pollutants[sampleTrain,]
dataTest <- pollutants[-sampleTrain,]

M0 <- lm(length ~ 1, data=pollutants)
Mfull <- lm(length ~ ., data=pollutants)
Mstart <- M0

MfwdAIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "forward", k = 2)
MfwdBIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "forward", k = log(N))
MbckAIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "backward", k = 2)
MbckBIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "backward", k = log(N))
MstepAIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "both", k = 2)
MstepBIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                trace = FALSE, direction = "both", k = log(N))

mspe <- function(model, df, y) {
  model.res <- y - predict(model, newdata=df)
  return(mean(model.res^2))
}

# cross validation
kfold <- 12
pollutants$index <- rep(1:kfold, each=N/kfold)

MfwdAIC.mspe <- rep(NA, kfold)
MfwdBIC.mspe <- rep(NA, kfold)
MbckAIC.mspe <- rep(NA, kfold)
MbckBIC.mspe <- rep(NA, kfold)
MstepAIC.mspe <- rep(NA, kfold)
MstepBIC.mspe <- rep(NA, kfold)

for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants$index != i)

  MfwdAIC.cv <- update(MfwdAIC, subset=train.ind)
  MfwdBIC.cv <- update(MfwdBIC, subset=train.ind)
  MbckAIC.cv <- update(MbckAIC, subset=train.ind)
  MbckBIC.cv <- update(MbckBIC, subset=train.ind)
  MstepAIC.cv <- update(MstepAIC, subset=train.ind)
  MstepBIC.cv <- update(MstepBIC, subset=train.ind)

  i_newdata <- pollutants[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  MfwdAIC.res <- mspe(MfwdAIC.cv, i_newdata, i_y)
  MfwdBIC.res <- mspe(MfwdBIC.cv, i_newdata, i_y)
  MbckAIC.res <- mspe(MbckAIC.cv, i_newdata, i_y)
  MbckBIC.res <- mspe(MbckBIC.cv, i_newdata, i_y)
  MstepAIC.res <- mspe(MstepAIC.cv, i_newdata, i_y)
  MstepBIC.res <- mspe(MstepBIC.cv, i_newdata, i_y)

  MfwdAIC.mspe[i] <- MfwdAIC.res
  MfwdBIC.mspe[i] <- MfwdBIC.res
  MbckAIC.mspe[i] <- MbckAIC.res
  MbckBIC.mspe[i] <- MbckBIC.res
  MstepAIC.mspe[i] <- MstepAIC.res
  MstepBIC.mspe[i] <- MstepBIC.res
  
}


#====
# test the three methods
pollutants_t <- pollutants_original[-1]

# forward selection using AIC
fwdAIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MfwdAIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
       trace = FALSE, direction = "forward", k = 2)

  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MfwdAIC.res <- mspe(i_MfwdAIC, i_newdata, i_y)

  fwdAIC.mspe[i] <- i_MfwdAIC.res
}

# backforward elimination using AIC
bckAIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MbckAIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "backward", k = 2)
  
  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MbckAIC.res <- mspe(i_MbckAIC, i_newdata, i_y)
  
  bckAIC.mspe[i] <- i_MbckAIC.res
}

# stepwise selection using AIC
stepAIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MstepAIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "both", k = 2)
  
  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MstepAIC.res <- mspe(i_MstepAIC, i_newdata, i_y)
  
  stepAIC.mspe[i] <- i_MstepAIC.res
}

# forward selection using BIC
fwdBIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MfwdBIC <- step(object = M0, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "forward", k = log(length(train.ind)))
  
  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MfwdBIC.res <- mspe(i_MfwdBIC, i_newdata, i_y)
  
  fwdBIC.mspe[i] <- i_MfwdBIC.res
}

# backforward elimination using AIC
bckBIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MbckBIC <- step(object = Mfull, scope = list(lower = M0, upper = Mfull), 
                    trace = FALSE, direction = "backward", k = log(length(train.ind)))
  
  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MbckBIC.res <- mspe(i_MbckBIC, i_newdata, i_y)
  
  bckBIC.mspe[i] <- i_MbckBIC.res
}

# stepwise selection using BIC
stepBIC.mspe <- rep(NA, kfold)
for (i in 1:kfold) {
  print(paste0("i = ", i))
  train.ind <- which(pollutants$index != i)
  
  M0 <- lm(length ~ 1, data=pollutants_t[train.ind, ])
  Mstart <- M0
  Mfull <- lm(length ~ ., data=pollutants_t[train.ind, ])
  i_MstepBIC <- step(object = Mstart, scope = list(lower = M0, upper = Mfull), 
                     trace = FALSE, direction = "both", k = log(length(train.ind)))
  i_newdata <- pollutants_t[-train.ind, ]
  i_y <- pollutants$length[-train.ind]
  i_MstepBIC.res <- mspe(i_MstepBIC, i_newdata, i_y)
  
  stepBIC.mspe[i] <- i_MstepBIC.res
}
print(paste("fwd AIC MSPE:", mean(fwdAIC.mspe)))
print(paste("fwd BIC MSPE:", mean(fwdBIC.mspe)))
print(paste("step AIC MSPE:", mean(stepAIC.mspe)))
print(paste("bck AIC MSPE:", mean(bckAIC.mspe)))
print(paste("bck BIC MSPE:", mean(bckBIC.mspe)))
print(paste("step BIC MSPE:", mean(stepBIC.mspe)))

print(mean(MfwdAIC.mspe))
print(mean(MfwdBIC.mspe))
print(mean(MstepAIC.mspe))
print(mean(MbckAIC.mspe))
print(mean(MbckBIC.mspe))
print(mean(MstepBIC.mspe))
