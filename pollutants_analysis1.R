library(corrplot)
data_dir <-  r"(C:\Users\Frank Shi\Documents\FrankS\Waterloo\pollutants)"
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

