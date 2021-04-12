set.seed(123) ## for reproducibility
library(glmnet)
library(corrplot)
data_dir <-  "/home/logan/Desktop/STAT 331/Final/pollutants"
setwd(data_dir)
pollutants_original <- read.csv(file='pollutants.csv')

#===============================================
# plot correlation bewtween all covariates
excl_cols <- c('X', 'length')
pollutants1 <- pollutants_original[!colnames(pollutants_original) %in% excl_cols]
pollutants_corr <- cor(pollutants1)
corrplot(pollutants_corr, order="hclust", method="color", type="full", 
         tl.col = "black")

# take out a few covariates based on the correlation heat map
excl_cols1 <- c('X', paste0('POP_PCB', c(1, 2, 4, 5, 7, 9)),
                paste0('POP_furan', c(1, 2)), 'POP_dioxin1',
                'lymphocyte_pct')
pollutants2 <- pollutants_original[!colnames(pollutants_original) %in% excl_cols1]
pollutants_corr2 <- cor(pollutants2)
corrplot(pollutants_corr2, order="hclust", method="color", type="full", 
         tl.col = "black")


#===============================================
# after the first round of filtering of covariates, deal with categorical covariates
pollutants <- pollutants2
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


#===============================================
# group different age groups
pollutants['agegroup'] <- 0
pollutants[pollutants['ageyrs'] >= 20 & pollutants['ageyrs'] <= 35, 'agegroup'] = 1
pollutants[pollutants['ageyrs'] >= 36 & pollutants['ageyrs'] <= 50, 'agegroup'] = 2
pollutants[pollutants['ageyrs'] >= 51, 'agegroup'] = 3
agegroup_list <- c('young', 'middleage', 'elder')
agegroup <- agegroup_list[pollutants$agegroup]
pollutants$agegroup <- factor(agegroup, levels=agegroup_list)

# lastly, take out ageyrs because we already have agegroup
pollutants <- pollutants[!colnames(pollutants) %in% c('ageyrs')]

# the most basic steps of data manipulation complete

# some ad-hoc covariates exclusions
# adhoc_excl <- c('edu_cat', 'race_cat')
# pollutants <- pollutants[!colnames(pollutants) %in% adhoc_excl]

#===============================================
# establishing benchmarks
model_full <- lm(length ~ ., data=pollutants)
model_0 <- lm(length ~ 1, data=pollutants)
model_start_aic <- lm(length ~ 1, data=pollutants)
# aic stepwise
system.time({
  Mstep_aic <- step(object = model_start_aic,
                scope = list(lower = model_0, upper = model_full),
                direction = "both", 
                trace = 1, 
                k = 2)
})

# bic stepwise
model_start_bic <- lm(length ~ 1, data=pollutants)
system.time({
  Mstep_bic <- step(object = model_start_bic,
                scope = list(lower = model_0, upper = model_full),
                direction = "both", 
                trace = 1, 
                k = log(nrow(pollutants)))
})

#  lasso
y <- pollutants[['length']]
X <- model.matrix(model_full)
model_lasso <- glmnet(x=X, y=y, alpha=1, nlambda=500)
cv_fit_lasso <- cv.glmnet(x=X, y=y, alpha=1, nlambda=500)
print(coef(model_lasso, s="lambda.min"))
print(coef(cv_fit_lasso, s="lambda.min"))

# build the models with the covariates chosen by these model selection algorithms
stepwise_covariates <- c('agegroup', 'male', 'ln_lbxcot', 'yrssmoke')
stepwise_dataset <- pollutants[colnames(pollutants) %in% c('length', stepwise_covariates)]
model_step <- lm(length ~ ., data=stepwise_dataset)
lasso_covariates <- c('POP_PCB8', 'POP_dioxin2', 'POP_dioxin3', 'POP_furan3', 'monocyte_pct', 
                      'BMI', 'edu_cat', 'race_cat', 'male', 'yrssmoke', 'ln_lbxcot', 'agegroup')
lasso_dataset <- pollutants[colnames(pollutants) %in% c('length', lasso_covariates)]
model_lasso <- lm(length ~ ., data=lasso_dataset)

#===============================================
# trying interaction terms
# on stepwise covariates
stepwise_age_interact <- lm(length ~ . + agegroup * ln_lbxcot + agegroup * yrssmoke, data=stepwise_dataset)
