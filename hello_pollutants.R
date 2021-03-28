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
covariates <- colnames(pollutants)[!(colnames(pollutants) %in% c('length', pcb_colnames, dioxin_colnames, furan_colnames))]

# fit the response variable to the pollutant exposures only, for reference
exposures_only <- pollutants[c('length', pcb_colnames, dioxin_colnames, furan_colnames)]
exposure_model <- lm(length ~ ., data=exposures_only)

# fit the pollutant exposures to the external covariates one by one
for (i in 1:pcb_kinds) {
  i_set <- pollutants[c(paste0('POP_PCB', i), covariates)]
  assign(paste0('pcb_', i), lm(as.formula(paste0('POP_PCB', i, '~.')), data=i_set))
}

# dioxin exposure
for (i in 1:dioxin_kinds) {
  i_set <- pollutants[c(paste0('POP_dioxin', i), covariates)]
  assign(paste0('dioxin_', i), lm(as.formula(paste0('POP_dioxin', i, '~.')), data=i_set))
}

# dioxin exposure
for (i in 1:furan_kinds) {
  i_set <- pollutants[c(paste0('POP_furan', i), covariates)]
  assign(paste0('furan_', i), lm(as.formula(paste0('POP_furan', i, '~.')), data=i_set))
}

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


