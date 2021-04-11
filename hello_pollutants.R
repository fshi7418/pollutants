data_dir <-  r"(C:\Users\Frank Shi\Documents\FrankS\Waterloo\pollutants)"
setwd(data_dir)
pollutants <- read.csv(file='pollutants.csv')

# some constants
pcb_kinds <- 11
dioxin_kinds <- 3
furan_kinds <- 4

# make a few things categorical
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

# get rid of the first column
pollutants <- pollutants[-1]

# fit the full model just for reference
full_model <- lm(length ~ ., data=pollutants)

# fit some models on only some columns
pcb_colnames <- paste('POP_PCB', 1:pcb_kinds, sep='')
dioxin_colnames <- paste('POP_dioxin', 1:dioxin_kinds, sep='')
furan_colnames <- paste('POP_furan', 1:furan_kinds, sep='')

# all the other column names. call them external covariates
covariates_ext <- colnames(pollutants)[!(colnames(pollutants) %in% c('length', pcb_colnames, dioxin_colnames, furan_colnames))]
# exclude a few covariates that intuitively might not contribute to the response
covariates_excl <- c('edu_cat', 'race_cat')
covariates_ext2 <- covariates_ext[!(covariates_ext %in% covariates_excl)]

# fit the response variable to the pollutant exposures only, for reference
exposures_only <- pollutants[c('length', pcb_colnames, dioxin_colnames, furan_colnames)]
exposure_model <- lm(length ~ ., data=exposures_only)


# relationship between the organic pollutants
# compute VIFs for each pollutant while regressing against others
# i.e. regressing PCB_1 against PCB2--11 and record VIF_1 for PCB_1, and so on;
#      regressiong dioxin_1 against dioxin2--3 and record VIF_1, for dioxin_1, and so on...
pcb_vifs <- c()
pcb_set <- pollutants[pcb_colnames]
for (i in 1:pcb_kinds) {
  i_model <- lm(as.formula(paste0('POP_PCB', i, '~.')), data=pcb_set)
  vif_i <- 1 / (1 - summary(i_model)$r.squared)
  pcb_vifs <- c(pcb_vifs, vif_i)
}

# for dioxins
dioxin_vifs <- c()
dioxin_set <- pollutants[dioxin_colnames]
for (i in 1:dioxin_kinds) {
  i_model <- lm(as.formula(paste0('POP_dioxin', i, '~.')), data=dioxin_set)
  vif_i <- 1 / (1 - summary(i_model)$r.squared)
  dioxin_vifs <- c(dioxin_vifs, vif_i)
}

# for furans
furan_vifs <- c()
furan_set <- pollutants[furan_colnames]
for (i in 1:furan_kinds) {
  i_model <- lm(as.formula(paste0('POP_furan', i, '~.')), data=furan_set)
  vif_i <- 1 / (1 - summary(i_model)$r.squared)
  furan_vifs <- c(furan_vifs, vif_i)
}


# one-by-one elimination of PCBs based on VIF
remaining_pcb_colnames <- pcb_colnames
pcb_set_iter <- pcb_set
for (i in 1:pcb_kinds) {  # the elimnation is performed at most p - 1 times
  num_covariates <- dim(pcb_set_iter)[2]
  delete_flag <- F
  for (j in 1:num_covariates) { # fit each covariates against the other covariates
    xj_name <- colnames(pcb_set_iter)[j]  # get the j-th covariate name
    j_m <- lm(as.formula(paste(xj_name, '~ .')), data=pcb_set_iter) # fit the model
    vif_j <- 1 / (1 - summary(j_m)$r.squared) # compute VIF_j
    print(paste(xj_name, 'VIF_j =', vif_j))
    if (vif_j > 10) {
      delete_flag <- T
      print(paste('deleting', xj_name))
      break
    }
  }
  if (delete_flag) {
    # take out the j-th covariate from the columns 
    pcb_set_iter <- pcb_set_iter[colnames(pcb_set_iter)[-j]]
    # print(paste('colnames af deletion:', colnames(covariates_training)))
  } else {
    # if all the VIFs are smaller than 10, return all the columns and break out of the loop
    chosen_pcbs <- colnames(pcb_set_iter)
    break
  }
}


pcb_set <- pcb_set[chosen_pcbs]
pcb_kinds <- dim(pcb_set)[2]

# fit the pollutant exposures to the external covariates one by one
for (pcb in chosen_pcbs) {
  i <- gsub('POP_PCB', '', pcb)
  i_set <- pollutants[c(pcb, covariates_ext2)]
  assign(paste0('pcb_', i), lm(as.formula(paste0('POP_PCB', i, '~.')), data=i_set))
}

# dioxin exposure
for (i in 1:dioxin_kinds) {
  i_set <- pollutants[c(paste0('POP_dioxin', i), covariates_ext2)]
  assign(paste0('dioxin_', i), lm(as.formula(paste0('POP_dioxin', i, '~.')), data=i_set))
}

# dioxin exposure
for (i in 1:furan_kinds) {
  i_set <- pollutants[c(paste0('POP_furan', i), covariates_ext2)]
  assign(paste0('furan_', i), lm(as.formula(paste0('POP_furan', i, '~.')), data=i_set))
}


# soe adhod investigations
yrssmoke_v_ageyrs <- lm(pollutants$yrssmoke ~ pollutants$ageyrs)
male_v_yrssmoke <- lm(pollutants$yrssmoke ~ pollutants$male)
yrssmoke_v_cotinine <- lm(pollutants$yrssmoke ~ pollutants$ln_lbxcot)


# run LASSO on the pollutant models
set.seed(123) ## for reproducibility
library(glmnet)

# check if the lasso variables for each pollutant agrees with the single-variate cases
for (pcb in chosen_pcbs) {
  print(paste0('current PCB is ', pcb))
  i <- gsub('POP_PCB', '', pcb)
  i_set <- pollutants[c(pcb, covariates_ext2)]
  i_lm <- lm(as.formula(paste0(pcb, '~.')), data=i_set)
  i_X <- model.matrix(i_lm)
  i_y <- pollutants[[pcb]]
  i_lasso <- glmnet(x=i_X, y=i_y, 
                    alpha=1,
                    nlambda=500)
  i_cvfit_lasso <- cv.glmnet(x=i_X, y=i_y, 
                             alpha=1,
                             nlambda=500)
  print(coef(i_cvfit_lasso, s = "lambda.min"))## alternatively could use "lambda.1se"
}

for (i in 1:furan_kinds) {
  i_name <- paste0('POP_furan', i)
  print(paste0('current furan is ', i_name))
  i_set <- pollutants[c(i_name, covariates_ext2)]
  i_lm <- lm(as.formula(paste0(i_name, '~.')), data=i_set)
  i_X <- model.matrix(i_lm)
  i_y <- pollutants[[i_name]]
  i_lasso <- glmnet(x=i_X, y=i_y, 
                    alpha=1,
                    nlambda=500)
  i_cvfit_lasso <- cv.glmnet(x=i_X, y=i_y, 
                             alpha=1,
                             nlambda=500)
  print(coef(i_cvfit_lasso, s = "lambda.min"))## alternatively could use "lambda.1se"
}

for (i in 1:dioxin_kinds) {
  i_name <- paste0('POP_dioxin', i)
  print(paste0('current dioxin is ', i_name))
  i_set <- pollutants[c(i_name, covariates_ext2)]
  i_lm <- lm(as.formula(paste0(i_name, '~.')), data=i_set)
  i_X <- model.matrix(i_lm)
  i_y <- pollutants[[i_name]]
  i_lasso <- glmnet(x=i_X, y=i_y, 
                    alpha=1,
                    nlambda=500)
  i_cvfit_lasso <- cv.glmnet(x=i_X, y=i_y, 
                             alpha=1,
                             nlambda=500)
  print(coef(i_cvfit_lasso, s = "lambda.min"))## alternatively could use "lambda.1se"
}