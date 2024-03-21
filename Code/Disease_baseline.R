# File to read icd code for spluiting in new columns: 

hesin_diag <- read.table(file = 'C:/Users/User/Desktop/Data/hesin_diag.txt' ,header = TRUE , sep = "\t" )%>% 
  select('requires.eid','ins_index','diag_icd10')

# Metformin supplement from Medication dataframe:
# medication file represents seperate file of Medications of Liver diseases used. 

# Separates rows with Metformin supplements : 
metformin_rows <- medication[rowSums(medication[, grepl("Medications_", names(medication))] == "metformin", na.rm = TRUE) > 0, ]

# metformin_df consists columns with Metfomin Eids 

metformin_df <- metformin_rows[, c('eid', grep("Medications_", names(medication), value = TRUE))]
metformin_df = metformin_df[!duplicated(metformin_df), ]

# data table of metformin_df :
metformin <- data.table(metformin_df)

# Separates rows with pioglitazone supplements : 
pioglitazone_rows <- medication[rowSums(medication[, grepl("Medications_", names(medication))] == "pioglitazone", na.rm = TRUE) > 0, ]

# metformin_df consists columns with pioglitazone Eids 

pioglitazone_df <- pioglitazone_rows[, c('eid', grep("Medications_", names(medication), value = TRUE))]
pioglitazone = pioglitazone_df[!duplicated(pioglitazone_df), ]

# Separates rows with ramipril supplements : 
ramipril_rows <- medication[rowSums(medication[, grepl("Medications_", names(medication))] == "ramipril", na.rm = TRUE) > 0, ]

# metformin_df consists columns with ramipril Eids 

ramipril_df <- ramipril_rows[, c('eid', grep("Medications_", names(medication), value = TRUE))]
ramipril_df = ramipril_df[!duplicated(ramipril_df), ]

# data table of ramipril_df :
ramipril <- data.table(ramipril_df)

# Identify rows that are not related to metformin and create negative control dataframe 
non_metformin <- medication[!medication$eid %in% metformin$eid, ]

# Create data table for negative control
non_metformin = non_metformin[!duplicated(non_metformin), ]
non_metformin <- data.table(non_metformin)


#MASLD icd code K760 to create new df. 

MASLD <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K760'))
MASLD <- as.data.frame(MASLD[!duplicated(MASLD$eid_1), ])
MASLD$Diagnosis <- 'K760'
MASLD <- select(MASLD, c(eid_1, Diagnosis))
setnames(MASLD,"eid_1","eid")

# NASH dataframe K758 ICD CODE: 

NASH <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K758'))
NASH <- as.data.frame(NASH[!duplicated(NASH$eid_1), ])
NASH$Diagnosis <- 'K758'
NASH <- select(NASH, c(eid_1, Diagnosis))
setnames(NASH,"eid_1","eid")

# hypertension :

hypertension_fr <- as.data.frame(hypertension_fr[!duplicated(hypertension_fr$eid_1), ]) # duplicates removal
hypertension_fr <- select(hypertension_fr, c(eid_1, diag_icd10)) # selecting Columns for hypertension
setnames(hypertension_fr,"eid_1","eid") # set column names
hypertension_fr$hypertension<-replace_na(hypertension_fr$hypertension, 0) # replace NAs of hypertension

# liver disease :
liverdisease <- as.data.frame(liverdisease[!duplicated(liverdisease$eid), ])
liverdisease <- setnames(liverdisease,"eid_1","eid")

# MASLD metformin percentage calculation against metformin df:

MASLD_metformin_percent <- length(intersect(metformin_MASLD$eid, metformin_df$eid)) / length(metformin_df$eid) * 100

# IDs to use

metformin$metformin <-1
MASLD$MASLD <-1
NASH$nash <-1
pioglitazone$pioglitazone <-1
ramipril$ramipril <-1
hypertension_fr$hypertension <-1
death_df$death <-1
liverdisease$liverdisease <- 1 

# ALL file merge with metformin for MASLD and NASH :

MASLD_merge <-merge(baseline_df, MASLD, by.x="eid", all.x = TRUE, by.y="eid", all.y = TRUE )
metformin_MASLD <-merge(MASLD_merge, metformin, by.x="eid", all.x = TRUE) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin', 'MASLD')
ALL <-merge(metformin_MASLD, NASH, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE )  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin', 'MASLD', 'nash')

# pioglitazone merge with MASLD and NASH :
pioglitazone_MASLD <-merge(MASLD_merge, pioglitazone, by.x="eid", all.x = TRUE) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','pioglitazone', 'MASLD')
ALL_pioglitazone <-merge(pioglitazone_MASLD, NASH, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','pioglitazone', 'MASLD', 'nash')

# Ramipril merge with ALL database :
ramipril_MASLD <-merge(MASLD_merge, ramipril, by.x="eid", all.x = TRUE) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','ramipril', 'MASLD')
ALL_ramipril <-merge(ramipril_MASLD, NASH, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','ramipril', 'MASLD', 'nash')
ALL_ramipril <-merge(ALL_ramipril, hypertension_fr, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE) %>%
  select('eid','Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','ramipril', 'MASLD', 'nash', 'hypertension')

# death df from MASLD and NASH death of petinsts : 
ALL_Death <-merge(ALL, death_df, by.x="eid", by.y="eid_1", all.x = TRUE, all.y = TRUE )  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin','death')
ALL_Death <-merge(ALL_Death, liverdisease, by.x="eid", by.y="eid", all.x = TRUE )  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin','death','liverdisease')

# Nutrition associated analysis :
ALL_Nutrition <- merge(nutritrion_df, liverdisease, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE)  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','liverdisease','meals', 'spz_diet', 'Alcohol_int_10y', 'Alcohol_int_fr','Alcohol_cons')

# Reshaping metformin from different value to 0:

ALL_ramipril$Diabetes <- ifelse(ALL$Diabetes == 1, 1, 0)

# Perform t-test for metformin patients: 

MASLD_ttest_metformin <- t.test(ALL$metformin ~ ALL$MASLD)
NASH_t_test_metformin <- t.test(ALL$metformin ~ ALL$nash)
Gender_t_test_metformin <- t.test(ALL$metformin ~ ALL$Gender)
diabetes_t_test_metformin <- t.test(ALL$metformin ~ ALL$Diabetes)
BMI_t_test_metformin <- t.test(ALL$metformin , ALL$BMI)
MASLD_t_test <- t.test(ALL$metformin, ALL$MASLD)
NASH_t_test <- t.test(ALL$metformin, ALL$nash)

# NAs to transform into 0 for pioglitazone medication:

pioglitazone_MASLD$MASLD<-replace_na(pioglitazone_MASLD$MASLD, 0)
pioglitazone_MASLD$pioglitazone<-replace_na(pioglitazone_MASLD$pioglitazone, 0)
ALL_pioglitazone$nash<-replace_na(ALL_pioglitazone$nash, 0)
ALL_pioglitazone$MASLD<-replace_na(ALL_pioglitazone$MASLD, 0)
ALL_pioglitazone$pioglitazone<-replace_na(ALL_pioglitazone$pioglitazone, 0)
ALL_pioglitazone$Diabetes[ALL_pioglitazone$Diabetes == -3 | ALL_pioglitazone$Diabetes == -1 | is.na(ALL_pioglitazone$Diabetes)] <- 0

#Pioglitazone t test for p value

MASLD_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$MASLD)
NASH_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$nash)
Gender_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$Gender)
diabetes_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone ,ALL_pioglitazone$Diabetes)
BMI_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$BMI)

# NAs to transform into 0 for ramipril medication:

ramipril_MASLD$MASLD<-replace_na(ramipril_MASLD$MASLD, 0)
ramipril_MASLD$ramipril<-replace_na(ramipril_MASLD$ramipril, 0)
ALL_ramipril$nash<-replace_na(ALL_ramipril$nash, 0)
ALL_ramipril$MASLD<-replace_na(ALL_ramipril$MASLD, 0)
ALL_ramipril$ramipril<-replace_na(ALL_ramipril$ramipril, 0)
ALL_ramipril$hypertension<-replace_na(ALL_ramipril$hypertension, 0)
ALL_ramipril$Diabetes[ALL_ramipril$Diabetes == -3 | ALL_ramipril$Diabetes == -1 | is.na(ALL_ramipril$Diabetes)] <- 0

# for metformin replace NAs to 0 :

ALL$MASLD<-replace_na(ALL$MASLD, 0)
ALL$metformin<-replace_na(ALL$metformin, 0)
ALL$nonmasld<-replace_na(ALL$nonmasld, 0)
ALL$nash<-replace_na(ALL$nash, 0)
ALL$Diabetes <-replace_na(ALL$Diabetes, 0)
ALL[is.na(ALL)] <- 0
ALL$Diabetes[ALL$Diabetes == -3 | ALL$Diabetes == -1 | is.na(ALL$Diabetes)] <- 0

# Nutrition NAs :
ALL_Nutrition[is.na(ALL_Nutrition)] <- 0

# death NAs to convert 
ALL_Death$death <-replace_na(ALL_Death$death, 0) # death df 
ALL_Death$liverdisease <-replace_na(ALL_Death$liverdisease, 0)

# Metformin group by 

metformin_group <- ALL[!is.na(ALL$metformin), ]
nonmetformin_group <- ALL[!is.na(ALL$nonmetformin), ]

# MASLD with Non metformin : 

common8 <- intersect(non_metformin$eid, MASLD$eid)
non_metformin_common <- non_metformin[non_metformin$eid %in% common8, ]
MASLD_common <- MASLD[MASLD$eid %in% common8, ]
non_metformin_common$eid <- as.integer(non_metformin_common$eid)
MASLD_non_metformin <- merge(non_metformin_common, MASLD_common, 
                        by.x = "eid", by.y = "eid", suffixes = c('_Non_metformin', '_MASLD'))
MASLD_non_metformin <- MASLD_non_metformin %>% 
  select(eid, Diagnosis)

# MASLD metformin percentage calculation against non-metformin df:

MASLD_nonmetformin_percent <- length(intersect(MASLD_non_metformin$eid, non_metformin$eid)) / length(non_metformin$eid) * 100

# Mean and standard deviation with/without metformin dataframes:

mean_MASLD_metformin <- mean(metformin_MASLD$eid, na.rm = TRUE)
mean_MASLD_non_metformin <- mean(MASLD_non_metformin$eid, na.rm = TRUE)
sd_MASLD_metformin<- sd(metformin_MASLD$eid, na.rm = TRUE)
sd_MASLD_non_metformin <- sd(MASLD_non_metformin$eid, na.rm = TRUE)

# adding dataframe numeric fomrmat:  
#mean for metformin df:
mean_MASLD_metformin <- as.numeric(mean_MASLD_metformin)
mean_MASLD_non_metformin <- as.numeric(mean_MASLD_non_metformin)
sd_MASLD_metformin <- as.numeric(sd_MASLD_metformin)
sd_MASLD_non_metformin <- as.numeric(sd_MASLD_non_metformin)

# standardized mean difference for MASLD ICD 10 with/without metformin:

smd_MASLD <- (mean_MASLD_metformin - mean_MASLD_non_metformin) / sqrt((sd_MASLD_metformin^2 + sd_MASLD_non_metformin^2) / 2)

# Ramipril :

matched31 <- matched3[matched3$ramipril == 1, ]
matched32 <- matched3[matched3$ramipril == 0, ]

# Calculate mean and standard deviation of age :
ramipril_age1 <- mean(matched31$Age_AC)
ramipril_age2 <- sd(matched31$Age_AC)
ramipril_age3 <- mean(matched32$Age_AC)
ramipril_age4 <- sd(matched32$Age_AC)

# Calculate mean and standard deviation of BMI : 
ramipril_BMI1 <- mean(matched31$BMI,na.rm = TRUE)
ramipril_BMI2 <- sd(matched31$BMI,na.rm = TRUE)
ramipril_BMI3 <- mean(matched32$BMI,na.rm = TRUE)
ramipril_BMI4 <- sd(matched32$BMI,na.rm = TRUE)

# Calculate standardized mean difference (SMD) for age
sd_age_ramipril <- sqrt((ramipril_age2^2 + ramipril_age4^2) / 2)
smd_age_ramipril <- (ramipril_age1 - ramipril_age3) / sd_age_ramipril

# Calculate standardized mean difference (SMD) for BMI
sd_bmi_ramipril <- sqrt((ramipril_BMI2^2 + ramipril_BMI4^2) / 2)
smd_bmi_ramipril <- (ramipril_BMI1 - ramipril_BMI3) / sd_bmi_ramipril

# Male female calculation ramipril intake:

ramipril_male <- sum(matched31$Gender == 1)
ramipril_female <- sum(matched31$Gender == 0)
ramipril_total_individuals <- ramipril_male + ramipril_female
ramipril_percentage_male <- (ramipril_male / ramipril_total_individuals) * 100
ramipril_percentage_female <- (ramipril_female / ramipril_total_individuals) * 100

# Male female calculation non ramipril intake:

nonramipril_male <- sum(matched32$Gender == 1)
nonramipril_female <- sum(matched32$Gender == 0)
nonramipril_individuals <- nonramipril_male + nonramipril_female
nonramipril_percen_male <- (nonramipril_male / nonramipril_individuals) * 100
nonramipril_percen_female <- (nonramipril_female / nonramipril_individuals) * 100

# sum for gender of SMD. 
sum_gender1 <- sum(matched3$Gender == 0, na.rm = TRUE)
sum_gender2 <- sum(matched3$Gender == 1, na.rm = TRUE)

# Calculate proportions
prop_gender1 <- sum_gender1 / length(matched3$Gender)
prop_gender2 <- sum_gender2 / length(matched3$Gender)

# Calculate sample sizes
n1 <- sum(!is.na(matched3$Gender) & matched3$Gender == 0)
n2 <- sum(!is.na(matched3$Gender) & matched3$Gender == 1)

n1_numeric <- as.numeric(n1)
n2_numeric <- as.numeric(n2)

# Calculate standardized mean difference (SMD) for gender
ramipril_smd_gender <- (prop_gender1 - prop_gender2) / sqrt((n1_numeric * n2_numeric) / (n1_numeric + n2_numeric))

# ramipril consuming percentage of individuals with and without diabetes :  

ramipril_diabetes_yes <- sum(matched31$Diabetes == 1, na.rm = TRUE)
ramipril_diabetes_no <- sum(matched31$Diabetes == 0, na.rm = TRUE)
non_ramipril_diabetes_yes <- sum(matched32$Diabetes == 1, na.rm = TRUE)
non_ramipril_diabetes_no <- sum(matched32$Diabetes == 0, na.rm = TRUE)

# total number of diabetes individuals

tl_ind <- ramipril_diabetes_yes + ramipril_diabetes_no
tl_ind2 <- non_ramipril_diabetes_yes + non_ramipril_diabetes_no 

# Calculate the percentage of individuals with and without diabetes:

ramipril_percent <- (ramipril_diabetes_yes / tl_ind) * 100
non_ramipril_percent <- (non_ramipril_diabetes_yes / tl_ind2) * 100

# Calculate the sum of binary values for Diabetes being 1
sum_diabetes1 <- sum(matched3$Diabetes == 1, na.rm = TRUE)

# Calculate the proportion of Diabetes being 1
prop_diabetes1 <- sum_diabetes1 / length(matched1$Diabetes)

# Calculate the sample size for Diabetes being 1
n1_diabetes <- sum(!is.na(matched1$Diabetes) & matched1$Diabetes == 1)

# Calculate standardized mean difference (SMD) for Diabetes being 1
ramipril_smd_diabetes1 <- prop_diabetes1

# ramipril consuming MASLD patients : 
ramipril_MASLD_yes <- sum(matched_ramipril$ramipril == 1 & matched_ramipril$MASLD == 1)
ramipril_MASLD_NO <- sum(matched_ramipril2$ramipril == 0 & matched_ramipril2$MASLD == 1)

# number of individuals consuming ramipril with MASLD
ramipril_individuals <- nrow(matched_ramipril)
ramipril_ind <- nrow(matched_ramipril2)

# Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
ramipril_MASLD_percent <- (ramipril_MASLD_yes / ramipril_individuals) * 100
non_ramipril_MASLD_percent <- (ramipril_MASLD_NO / ramipril_ind) * 100

# Calculate the sum of binary values for Diabetes being 1
ramipril_MASLD_sum <- sum(matched3$MASLD == 1, na.rm = TRUE)

# Calculate the proportion of Diabetes being 1
ramipril_MASLD_prop <- ramipril_MASLD_sum / length(matched_ramipril$MASLD)

# Calculate the sample size for Diabetes being 1
ramipril_MASLD_n1 <- sum(!is.na(matched1$MASLD) & matched1$MASLD == 1)

# Calculate standardized mean difference (SMD) for Diabetes being 1
ramipril_smd_MASLD <- ramipril_MASLD_prop

# ramipril consuming MASLD patients : 
ramipril_MASH_yes <- sum(matched_ramipril$ramipril == 1 & matched_ramipril$nash == 1)
ramipril_MASH_NO <- sum(matched_ramipril2$ramipril == 0 & matched_ramipril2$nash == 1)

# Calculate the percentage of individuals where both ramipril and MASH are equal to 1
ramipril_MASH_percent <- (ramipril_MASH_yes / ramipril_individuals) * 100
non_ramipril_MASH_percent <- (ramipril_MASH_NO / ramipril_ind) * 100

#  the sum of binary values for MASH being 1
ramipril_MASH_sum <- sum(matched3$nash == 1, na.rm = TRUE)

# Calculate the proportion of MASH being 1
ramipril_MASH_smd <- ramipril_MASH_sum / length(matched1$nash)

# t test for ramipril

ramipril_t_test_MASLD <- t.test(matched3$ramipril , matched3$MASLD, na.rm = TRUE)
ramipril_t_test_MASH <- t.test(matched3$ramipril , matched3$nash, na.rm = TRUE)
ramipril_t_test_Diabetes <- t.test(matched3$Diabetes, matched3$ramipril , na.rm = TRUE)
ramipril_t_test_Gender <- t.test(matched3$ramipril , matched3$Gender, na.rm = TRUE)


# pioglitazone  :

matched21 <- matched2[matched2$pioglitazone == 1, ]
matched22 <- matched2[matched2$pioglitazone == 0, ]

# Calculate mean and standard deviation of age :
pioglitazone_age1 <- mean(matched21$Age_AC)
pioglitazone_age2 <- sd(matched21$Age_AC)
pioglitazone_age3 <- mean(matched22$Age_AC)
pioglitazone_age4 <- sd(matched22$Age_AC)

# Calculate mean and standard deviation of BMI : 
pioglitazone_BMI1 <- mean(matched21$BMI,na.rm = TRUE)
pioglitazone_BMI2 <- sd(matched21$BMI,na.rm = TRUE)
pioglitazone_BMI3 <- mean(matched22$BMI,na.rm = TRUE)
pioglitazone_BMI4 <- sd(matched22$BMI,na.rm = TRUE)

# Pooled standard deviation for age
pooled_sd_age <- sqrt((pioglitazone_age2^2 + pioglitazone_age4^2) / 2)
pioglitazone_smd_age <- (pioglitazone_age1 - pioglitazone_age3) / pooled_sd_age #Standardized mean difference for age

# Pooled standard deviation for BMI
pooled_sd_bmi <- sqrt((pioglitazone_BMI2^2 + pioglitazone_BMI4^2) / 2)
pioglitazone_smd_bmi <- (pioglitazone_BMI1 - pioglitazone_BMI3) / pooled_sd_bmi # Standardized mean difference for BMI

# Male female calculation pioglitazone intake:

num_male <- sum(matched21$Gender == 1)
num_female <- sum(matched21$Gender == 0)
total_individuals <- num_male + num_female
pioglitazone_percentage_male <- (num_male / total_individuals) * 100
pioglitazone_percentage_female <- (num_female / total_individuals) * 100

# Male female calculation non pioglitazone intake:

num_male <- sum(matched22$Gender == 1)
num_female <- sum(matched22$Gender == 0)
total_individuals <- num_male + num_female
pioglitazone_percentage_male <- (num_male / total_individuals) * 100
pioglitazone_percentage_female <- (num_female / total_individuals) * 100

# sum for gender of SMD. 
sum_gender1 <- sum(matched2$Gender == 0, na.rm = TRUE)
sum_gender2 <- sum(matched2$Gender == 1, na.rm = TRUE)

# Calculate proportions
prop_gender1 <- sum_gender1 / length(matched2$Gender)
prop_gender2 <- sum_gender2 / length(matched2$Gender)

# Calculate sample sizes
n1 <- sum(!is.na(matched2$Gender) & matched2$Gender == 0)
n2 <- sum(!is.na(matched2$Gender) & matched2$Gender == 1)

n1_numeric <- as.numeric(n1)
n2_numeric <- as.numeric(n2)

# Calculate standardized mean difference (SMD) for gender
ramipril_smd_gender <- (prop_gender1 - prop_gender2) / sqrt((n1_numeric * n2_numeric) / (n1_numeric + n2_numeric))

# pioglitazone consuming percentage of individuals with and without diabetes :  

pioglitazone_diabetes_yes <- sum(matched2$Diabetes == 1, na.rm = TRUE)
pioglitazone_diabetes_no <- sum(matched2$Diabetes == 0, na.rm = TRUE)
non_pioglitazone_diabetes_yes <- sum(matched22$Diabetes == 1, na.rm = TRUE)
non_pioglitazone_diabetes_no <- sum(matched22$Diabetes == 0, na.rm = TRUE)

# Calculate the sum of binary values for Diabetes being 1
pioglitazone_sum_diabetes1 <- sum(matched2$Diabetes, na.rm = TRUE)

# Calculate the proportion of Diabetes being 1
pioglitazone_diabetes1 <- pioglitazone_diabetes_yes / length(matched2$Diabetes)

# Calculate the sample size for Diabetes being 1
pioglitazone_n1_diabetes <- sum(!is.na(matched2$Diabetes) & matched2$Diabetes == 1)

# Calculate standardized mean difference (SMD) for Diabetes being 1
pioglitazone_smd_diabetes1 <- pioglitazone_diabetes1

# total number of diabetes individuals

tl_ind <- pioglitazone_diabetes_yes + pioglitazone_diabetes_no
tl_ind2 <- non_pioglitazone_diabetes_yes + non_pioglitazone_diabetes_no 

# Calculate the percentage of individuals with and without diabetes:

pioglitazone_percent <- (pioglitazone_diabetes_yes / tl_ind) * 100
non_pioglitazone_percent <- (non_pioglitazone_diabetes_yes / tl_ind2) * 100

# pioglitazone consuming MASLD patients : 
pioglitazone_MASLD_yes <- sum(matched2$pioglitazone == 1 & matched21$MASLD == 1)
pioglitazone_MASLD_NO <- sum(matched22$pioglitazone == 0 & matched22$MASLD == 1)

# number of individuals consuming pioglitazone with MASLD
pioglitazone_individuals <- nrow(matched21)
pioglitazone_ind <- nrow(matched22)

# Calculate the percentage of individuals where both pioglitazone and MASLD are equal to 1
pioglitazone_MASLD_percent <- (pioglitazone_MASLD_yes / pioglitazone_individuals) * 100
non_pioglitazone_MASLD_percent <- (pioglitazone_MASLD_NO / pioglitazone_ind) * 100

#  the sum of binary values for MASLD being 1
pioglitazone_MASLD_sum <- sum(matched2$MASLD == 1, na.rm = TRUE)
pioglitazone_MASLD_smd <- pioglitazone_MASLD_sum / length(matched2$MASLD)

# pioglitazone consuming MASLD patients : 
pioglitazone_MASH_yes <- sum(matched21$pioglitazone == 1 & matched21$nash == 1)
pioglitazone_MASH_NO <- sum(matched22$pioglitazone == 0 & matched22$nash == 1)

# Calculate the percentage of individuals where both pioglitazone and MASH are equal to 1
pioglitazone_MASH_percent <- (pioglitazone_MASH_yes / pioglitazone_individuals) * 100
non_pioglitazone_MASH_percent <- (pioglitazone_MASH_NO / pioglitazone_ind) * 100

#  the sum of binary values for MASH being 1
pioglitazone_MASH_sum <- sum(matched2$nash == 1, na.rm = TRUE)
pioglitazone_MASH_smd <- pioglitazone_MASH_sum / length(matched2$nash)

# t test for pioglitazone

pioglitazone_t_test_MASLD <- t.test(matched2$pioglitazone , matched2$MASLD, na.rm = TRUE)
pioglitazone_t_test_MASH <- t.test(matched2$pioglitazone , matched2$nash, na.rm = TRUE)
pioglitazone_t_test_Diabetes <- t.test(matched2$pioglitazone , matched2$Diabetes, na.rm = TRUE)
pioglitazone_t_test_Gender <- t.test(matched2$pioglitazone , matched2$Gender, na.rm = TRUE)

# T-Test at Age baseline T-Test with/without metformin:
#calculating Mean and SD of Metformin with Age of baseline :

# Calculate mean and standard deviation of BMI
matched1_metformin <- subset(matched1, metformin == 1)
matched1_nonmetformin <- subset(matched1, metformin == 0)

# BMI matched mean and SD for metformin and nonmetformin :
metformin_mean_bmi <- mean(matched1_metformin$BMI)
metformin_sd_bmi <- sd(matched1_metformin$BMI)
nonmetformin_mean_bmi <- mean(matched1_nonmetformin$BMI)
nonmetformin_sd_bmi <- sd(matched1_nonmetformin$BMI)

# Calculate pooled standard deviation
metformin_sd_bmi <- sqrt((metformin_sd_bmi^2 + nonmetformin_sd_bmi^2) / 2)
smd_bmi <- (metformin_mean_bmi - nonmetformin_mean_bmi) / metformin_sd_bmi


# Age matched mean and SD :

metformin_mean_age <- mean(matched1_metformin$Age_AC)
matformin_sd_age <- sd(matched1_metformin$Age_AC)
nonmetformin_mean_age <- mean(matched1_nonmetformin$Age_AC)
nonmetformin_sd_age <- sd(matched1_nonmetformin$Age_AC)

# BMI :
nonmetformin_sd_age <- sqrt((matformin_sd_age^2 + nonmetformin_sd_age^2) / 2)
nonmetformin_smd_age <- (matformin_sd_age - nonmetformin_sd_age) / nonmetformin_sd_age

# Create a subset where metformin is 1 and nash is 1
matched1_metformin_nash <- subset(matched1_metformin, nash == 1)
subset_count <- nrow(matched1_metformin_nash)
total_count <- nrow(matched1_metformin)
percentage <- (subset_count / total_count) * 100

# Create a subset where metformin is 1 and nash is 1
matched1_nonmetformin_nash <- subset(matched1_nonmetformin, nash == 1)
subset_count <- nrow(matched1_nonmetformin_nash)
total_count <- nrow(matched1_nonmetformin)
# Calculate the percentage
percentage <- (subset_count / total_count) * 100

# Create a subset where metformin is 1 and MASLD is 1
matched1_metformin_MASLD <- subset(matched1_metformin, MASLD == 1)
subset_count <- nrow(matched1_metformin_MASLD)
# Calculate the total number of rows in matched1
total_count <- nrow(matched1_metformin)
percentage <- (subset_count / total_count) * 100

# Create a subset where nonmetformin metformin and MASLD is 1
matched1_nonmetformin_MASLD <- subset(matched1_nonmetformin, MASLD == 1)
subset_count <- nrow(matched1_nonmetformin_MASLD)
# Calculate the total number of rows in matched1
total_count <- nrow(matched1_nonmetformin)
percentage <- (subset_count / total_count) * 100

# matched Diabetes
matched1_metformin_dianbetes <- subset(matched1_metformin)
subset_count <- nrow(matched1_metformin_dianbetes)
total_count <- nrow(matched1_nonmetformin)
# Calculate the percentage
percentage <- (subset_count / total_count) * 100

# Convert sample sizes to numeric
matched1_n1_diabetes_numeric <- as.numeric(matched1_n1_diabetes)
matched1_n2_diabetes_numeric <- as.numeric(matched1_n2_diabetes)

# Calculate standardized mean difference (SMD) for Diabetes
matched_smd_diabetes <- (matched_prop_diabetes1 - matched_prop_diabetes2) / sqrt((matched1_n1_diabetes_numeric * matched1_n2_diabetes_numeric) / (matched1_n1_diabetes_numeric + matched1_n2_diabetes_numeric))

# matched male metformin
matched1_male <- subset(matched1_metformin, Gender == 1 & metformin == 1)
subset_count2 <- nrow(matched1_male)
# Calculate the total number of rows in matched1
total_count <- nrow(matched1_metformin)
percentage <- (subset_count2 / total_count) * 100

# matched female metformin
matched1_female <- subset(matched1_metformin, Gender == 0 & metformin == 1)
subset_count2 <- nrow(matched1_female)
# Calculate the total number of rows in matched1
total_count <- nrow(matched1_metformin)
percentage <- (subset_count2 / total_count) * 100

# matched male non metformin
matched1_male <- subset(matched1_nonmetformin, Gender == 1 & metformin == 0)
subset_count2 <- nrow(matched1_male)
# Calculate the total number of rows in matched1
total_count <- nrow(matched1_metformin)
percentage <- (subset_count2 / total_count) * 100

# matched female 
matched1_female <- subset(matched1_nonmetformin, Gender == 0 & metformin == 0)
subset_count <- nrow(matched1_female)
# Calculate the total number of rows in matched1
total_count2 <- nrow(matched1_nonmetformin)
percentage <- (subset_count / total_count) * 100

# sum for gender of SMD. 
matched1_gender1 <- sum(matched1$Gender == 0, na.rm = TRUE)
matched1_gender2 <- sum(matched1$Gender == 1, na.rm = TRUE)

# Calculate proportions
matched_prop_gender1 <- matched1_gender1 / length(matched1$Gender)
matched_prop_gender2 <- matched1_gender2 / length(matched1$Gender)

# Calculate sample sizes
matched1_n1 <- sum(!is.na(matched1$Gender) & matched1$Gender == 0)
matched1_n2 <- sum(!is.na(matched1$Gender) & matched1$Gender == 1)

matched1_n1_numeric <- as.numeric(matched1_n1)
matched1_n2_numeric <- as.numeric(matched1_n2)

# Calculate standardized mean difference (SMD) for gender
matched_smd_gender <- (matched_prop_gender1 - matched_prop_gender2) / sqrt((matched1_n1_numeric * matched1_n2_numeric) / (matched1_n1_numeric + matched1_n2_numeric))

# sum for gender of SMD. 
mean_bmi <- mean(matched1$BMI, na.rm = TRUE)

# Calculate standard deviation for entire dataset
sd_bmi <- sd(matched1$BMI)

# Calculate Cohen's d
cohens_d <- mean_bmi / sd_bmi

# Calculate proportions MASLD
matched_prop_mafl1 <- sum(matched1$MASLD == 0, na.rm = TRUE) / length(matched1$MASLD)
matched_prop_mafl2 <- sum(matched1$MASLD == 1, na.rm = TRUE) / length(matched1$MASLD)

# Calculate sample sizes
matched1_n1_mafl <- sum(!is.na(matched1$MASLD) & matched1$MASLD == 0)
matched1_n2_mafl <- sum(!is.na(matched1$MASLD) & matched1$MASLD == 1)

# Convert sample sizes to numeric
matched1_n1_mafl_numeric <- as.numeric(matched1_n1_mafl)
matched1_n2_mafl_numeric <- as.numeric(matched1_n2_mafl)

# Calculate standardized mean difference (SMD) for MASLD
matched_smd_mafl <- (matched_prop_mafl1 - matched_prop_mafl2) / 
  sqrt((matched1_n1_mafl_numeric * matched1_n2_mafl_numeric) / 
         (matched1_n1_mafl_numeric + matched1_n2_mafl_numeric))

# Calculate proportions NASH
matched_prop_NASH1 <- sum(matched1$nash == 0, na.rm = TRUE) / length(matched1$nash)
matched_prop_NASH2 <- sum(matched1$nash == 1, na.rm = TRUE) / length(matched1$nash)

# Calculate sample sizes
matched1_n1_NASH <- sum(!is.na(matched1$nash) & matched1$nash == 0)
matched1_n2_NASH <- sum(!is.na(matched1$nash) & matched1$nash == 1)

# Convert sample sizes to numeric
matched1_n1_NASH_numeric <- as.numeric(matched1_n1_NASH)
matched1_n2_NASH_numeric <- as.numeric(matched1_n2_NASH)

# Calculate standardized mean difference (SMD) for NASH
matched_smd_NASH <- (matched_prop_NASH1 - matched_prop_NASH2) / 
  sqrt((matched1_n1_NASH_numeric * matched1_n2_NASH_numeric) / 
         (matched1_n1_NASH_numeric + matched1_n2_NASH_numeric))

# t test for metformin

matched1_t_test_MASLD <- t.test(matched1$metformin , matched1$MASLD, na.rm = TRUE)
matched1_t_test_MASH <- t.test(matched1$metformin , matched1$nash, na.rm = TRUE)
matched1_t_test_Diabetes <- t.test(matched1$metformin , matched1$Diabetes, na.rm = TRUE)
matched1_t_test_Gender <- t.test(matched1$metformin , matched1$Gender, na.rm = TRUE)
