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
MASLD <- as.data.frame(MASLD[!duplicated(MASLD$eid), ])
MASLD$Diagnosis <- 'K760'
MASLD <- select(MASLD, c(eid, Diagnosis))

# MASH dataframe K758 ICD CODE: 

MASH <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K758'))
MASH <- as.data.frame(MASH[!duplicated(MASH$eid), ])
MASH$Diagnosis <- 'K758'
MASH <- select(MASH, c(eid, Diagnosis))

# liver related Mortality

mortality_df <- death_df[death_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0'), ]

# Toxic liver:

Alcoholic_liver <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K70.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K70.0') %>% select(eid, icd_group)
Toxic_liver <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K71.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K71.0') %>% select(eid, icd_group)
Liver_failure <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K72.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K72.0') %>% select(eid, icd_group)
Liver_failure_unspecified <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K729')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K729') %>% select(eid, diag_icd10)
Chron_Hepatitis <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K73.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K73.0') %>% select(eid, icd_group)
fibrose_cirrho <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K74.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K74.0') %>% select(eid, icd_group)
Liverfibrose <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K740')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K740') %>% select(eid, diag_icd10)
biliary_cholangitis1 <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K743')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K743') %>% select(eid, diag_icd10)
biliary_cirrho2 <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K744')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K744') %>% select(eid, diag_icd10)
Biliary_cirrho <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K745')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K745') %>% select(eid, diag_icd10)
unspec_cirrho <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K746')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K746') %>% select(eid, diag_icd10)
inflamm_liver <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K75.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K75.0') %>% select(eid, icd_group)
Liver_abscess <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K750')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K750') %>% select(eid, diag_icd10)
otherD_liver <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K76.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K76.0') %>% select(eid, icd_group)
portal_hyper <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K766')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K766') %>% select(eid, diag_icd10)
LiverD_unspe <- hesin_diag %>% subset(startsWith(as.character(diag_icd10), 'K769')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(diag_icd10 = 'K769') %>% select(eid, diag_icd10)
Unclassified_Liverd <- hesin_diag %>% subset(startsWith(as.character(icd_group), 'K77.0')) %>% distinct(eid, .keep_all = TRUE) %>% mutate(icd_group = 'K77.0') %>% select(eid, icd_group)

# ALT covariates:

covariates <- covariates %>% select('eid','ALT','AST','GGT')
glucose <- read.csv(file = 'C:/Users/User/Desktop/Data/Labs.csv') %>% select('eid','X3Glucose')
setnames(glucose,"X3Glucose","GLS")
metacov <- merge(covariates, glucose, by.x="eid", all.x = TRUE) %>% select('eid','ALT','GLS','AST','GGT')

#ALL Liver disease:

All_liver <- bind_rows(Alcoholic_liver, Toxic_liver, Liver_failure, Chron_Hepatitis, fibrose_cirrho, inflamm_liver, otherD_liver, Unclassified_Liverd, MASLD, MASH)
All_liver <- All_liver %>% mutate(source = case_when(
    icd_group == 'K70.0' ~ 'Alcoholic_liver', icd_group == 'K71.0' ~ 'Toxic_liver', icd_group == 'K72.0' ~ 'Liver_failure',
    icd_group == 'K73.0' ~ 'Chron_Hepatitis', icd_group == 'K74.0' ~ 'fibrose_cirrho', icd_group == 'K75.0' ~ 'inflamm_liver',
    icd_group == 'K76.0' ~ 'otherD_liver', icd_group == 'K77.0' ~ 'Unclassified_Liverd', Diagnosis == 'K760' ~ 'MASLD', Diagnosis == 'K758' ~ 'MASH'))

All_liver[is.na(All_liver)] <- 0
All_liver <- as.data.frame(All_liver[!duplicated(All_liver$eid), ])
All_liver <- All_liver %>% select('eid','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','MASLD','MASH')

# hypertension :

hypertension_fr <- select(hypertension_fr, c(eid, diag_icd10)) # selecting Columns for hypertension
setnames(hypertension_fr,"eid_1","eid") # set column names
hypertension_fr$hypertension<-replace_na(hypertension_fr$hypertension, 0) # replace NAs of hypertension

# liver disease :

liverdisease <- as.data.frame(liverdisease[!duplicated(liverdisease$eid), ])
liverdisease <- liverdisease %>% select('eid','diag_icd10')

# IDs to use

metformin$metformin <-1
MASLD$MASLD <-1
MASH$MASH <-1
pioglitazone$pioglitazone <-1
ramipril$ramipril <-1
hypertension_fr$hypertension <-1
mortality_df$death <-1
liverdisease$liverdisease <- 1 
Alcoholic_liver$alcoholicliver <- 1
Toxic_liver$toxicliver <- 1
Liver_failure$Liverfailure <- 1
Chron_Hepatitis$ChronHepatitis <- 1
fibrose_cirrho$fibrosecirrho <- 1
inflamm_liver$inflammliver <- 1
otherD_liver$otherliverD <- 1
Unclassified_Liverd$UnclassLiverd <- 1

# ALL file merge with metformin for MASLD and MASH :

MASLD_merge <-merge(baseline_df, All_liver, by.x="eid", all.x = TRUE, by.y="eid", all.y = TRUE ) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')
ALL <-merge(MASLD_merge, metformin, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')
ALL <-merge(ALL, metacov, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT')
ALL <-merge(ALL, liverdisease, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease')
ALL <-merge(ALL, mortality_df, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
            'Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death')
ALL <-merge(ALL, mri, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
            'Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death','low_fat','high_fat')


# pioglitazone merge with MASLD and MASH :

ALL_pioglitazone <-merge(baseline_df, pioglitazone, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','pioglitazone')
ALL_pioglitazone <-merge(ALL_pioglitazone, All_liver, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE ) %>% select('eid','Age_AC', 'Gender','BMI','Smoking', 
         'Drinking', 'Qualifications', 'Diabetes','pioglitazone', 'MASLD', 'MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')
ALL_pioglitazone <-merge(ALL_pioglitazone, metacov, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
            'Diabetes','pioglitazone', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT')
ALL_pioglitazone <-merge(ALL_pioglitazone, liverdisease, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','pioglitazone', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease')
ALL_pioglitazone <-merge(ALL_pioglitazone, mortality_df, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','pioglitazone', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death')
ALL_pioglitazone <-merge(ALL_pioglitazone, mri, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
         'Diabetes','pioglitazone', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death','low_fat','high_fat')

# Ramipril merge with ALL database :

ramipril_baseline <-merge(baseline_df, ramipril, by.x="eid", all.x = TRUE) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','ramipril')
ALL_ramipril <-merge(All_liver, ramipril_baseline, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking',
       'Drinking', 'Qualifications', 'Diabetes','ramipril','MASLD', 'MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')
ALL_ramipril <-merge(ALL_ramipril, hypertension_fr, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking',
       'Drinking', 'Qualifications','hypertension', 'Diabetes','ramipril','MASLD', 'MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','hypertension')
ALL_ramipril <-merge(ALL_ramipril, metacov, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
       'hypertension','Diabetes','ramipril', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT')
ALL_ramipril <-merge(ALL_ramipril, liverdisease, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
      'hypertension','Diabetes','ramipril', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease')
ALL_ramipril <-merge(ALL_ramipril, mortality_df, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
      'Diabetes','ramipril','hypertension', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death')
ALL_ramipril <-merge(ALL_ramipril, mri, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
      'Diabetes','ramipril','hypertension', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','death','low_fat','high_fat')


# death df from MASLD and MASH death of patients : 

ALL_Death <-merge(ALL, death_df, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE )  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin','death')
ALL_Death <-merge(ALL_Death, ALL, by.x="eid", by.y="eid", all.x = TRUE )  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','metformin','death','liverdisease')

# ALT, AST, GLS at baseline 

ALL_normal <-merge(baseline_df, metformin, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
           'Diabetes','metformin')
ALL_normal <-merge(ALL_normal, metacov, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
           'Diabetes','metformin','ALT','GLS','AST','GGT')
ALL_normal <-merge(ALL_normal,All_liver, by.x="eid", all.x = TRUE) %>% select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications',
                                                                              'Diabetes','metformin','ALT','GLS','AST','GGT','MASLD','MASH')
ALL_normal$ALT <-replace_na(ALL_normal$ALT, 0)
ALL_normal$GLS <-replace_na(ALL_normal$GLS, 0)
ALL_normal$GGT <-replace_na(ALL_normal$GGT, 0)
ALL_normal$AST <-replace_na(ALL_normal$AST, 0)
ALL_normal$MASLD <-replace_na(ALL_normal$MASLD, 0)
ALL_normal$MASH <-replace_na(ALL_normal$MASH, 0)
ALL_normal$metformin <-replace_na(ALL_normal$metformin, 0)

ALL_normal <- ALL_normal %>% select('eid', 'Age_AC', 'Gender','BMI','Diabetes','metformin','ALT','GLS','AST','GGT','MASLD','MASH')
psm_model5 <- matchit(metformin ~ Age_AC + Gender + BMI + ALT + GLS + AST + GGT, data = ALL_normal, method = "nearest", ratio = 5)
summary(psm_model5)
matched5 <- match.data(psm_model5)

# Nutrition associated analysis :

ALL_Nutrition <- merge(nutritrion_df, liverdisease, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE) %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','liverdisease','meals', 'spz_diet', 'Alcohol_int_10y', 'Alcohol_int_fr','Alcohol_cons')


# NAs to transform into 0 for ramipril medication:

ALL_ramipril[is.na(ALL_ramipril)] <- 0
ALL_ramipril$Diabetes[ALL_ramipril$Diabetes == -3 | ALL_ramipril$Diabetes == -1 | is.na(ALL_ramipril$Diabetes)] <- 0
ALL_ramipril$ramipril <-replace_na(ALL_ramipril$ramipril, 0)

# for metformin replace NAs to 0 :

ALL[is.na(ALL)] <- 0
ALL$Diabetes[ALL$Diabetes == -3 | ALL$Diabetes == -1 | is.na(ALL$Diabetes)] <- 0
ALL$metformin <-replace_na(ALL$metformin, 0)
ALL$death <-replace_na(ALL$death, 0)

# Nutrition NAs :
ALL_Nutrition[is.na(ALL_Nutrition)] <- 0

# death NAs to convert 
ALL_Death$death <-replace_na(ALL_Death$death, 0) # death df 
ALL_Death$liverdisease <-replace_na(ALL_Death$liverdisease, 0)

# NAs to transform into 0 for pioglitazone medication:

ALL_pioglitazone[is.na(ALL_pioglitazone)] <- 0
ALL_pioglitazone$Diabetes[ALL_pioglitazone$Diabetes == -3 | ALL_pioglitazone$Diabetes == -1 | is.na(ALL_pioglitazone$Diabetes)] <- 0
ALL_pioglitazone$pioglitazone <-replace_na(ALL_pioglitazone$pioglitazone, 0)

# matched all medications, survival curve for liver related mortality
# AUC RUC curve kaplan maier 

# Perform t-test for metformin patients: 

MASLD_ttest_metformin <- t.test(ALL$metformin ~ ALL$MASLD)
MASH_t_test_metformin <- t.test(ALL$metformin ~ ALL$MASH)
Gender_t_test_metformin <- t.test(ALL$metformin ~ ALL$Gender)
diabetes_t_test_metformin <- t.test(ALL$metformin ~ ALL$Diabetes)
BMI_t_test_metformin <- t.test(ALL$metformin , ALL$BMI)
MASLD_t_test <- t.test(ALL$metformin, ALL$MASLD)
MASH_t_test <- t.test(ALL$metformin, ALL$MASH)

#Pioglitazone t test for p value

MASLD_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$MASLD)
MASH_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$MASH)
Gender_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$Gender)
diabetes_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone ,ALL_pioglitazone$Diabetes)
BMI_ttest_Pioglitazone <- t.test(ALL_pioglitazone$pioglitazone , ALL_pioglitazone$BMI)

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

# Calculate mean and standard deviation of ALT : 
ramipril_ALT1 <- mean(matched31$ALT,na.rm = TRUE)
ramipril_ALT2 <- sd(matched31$ALT,na.rm = TRUE)
ramipril_ALT3 <- mean(matched32$ALT,na.rm = TRUE)
ramipril_ALT4 <- sd(matched32$ALT,na.rm = TRUE)

# Calculate mean and standard deviation of GLS : 
ramipril_GLS1 <- mean(matched31$GLS,na.rm = TRUE)
ramipril_GLS2 <- sd(matched31$GLS,na.rm = TRUE)
ramipril_GLS3 <- mean(matched32$GLS,na.rm = TRUE)
ramipril_GLS4 <- sd(matched32$GLS,na.rm = TRUE)

# Calculate mean and standard deviation of AST : 
ramipril_AST1 <- mean(matched31$AST,na.rm = TRUE)
ramipril_AST2 <- sd(matched31$AST,na.rm = TRUE)
ramipril_AST3 <- mean(matched32$AST,na.rm = TRUE)
ramipril_AST4 <- sd(matched32$AST,na.rm = TRUE)

# Calculate mean and standard deviation of GGT : 
ramipril_GGT1 <- mean(matched31$GGT,na.rm = TRUE)
ramipril_GGT2 <- sd(matched31$GGT,na.rm = TRUE)
ramipril_GGT3 <- mean(matched32$GGT,na.rm = TRUE)
ramipril_GGT4 <- sd(matched32$GGT,na.rm = TRUE)

# Calculate mean and standard deviation of LFAT : 
ramipril_LFAT1 <- mean(matched31$low_fat,na.rm = TRUE)
ramipril_LFAT2 <- sd(matched31$low_fat,na.rm = TRUE)
ramipril_LFAT3 <- mean(matched32$low_fat,na.rm = TRUE)
ramipril_LFAT4 <- sd(matched32$low_fat,na.rm = TRUE)

# Calculate mean and standard deviation of MFAT : 
ramipril_MFAT1 <- mean(matched31$high_fat,na.rm = TRUE)
ramipril_MFAT2 <- sd(matched31$high_fat,na.rm = TRUE)
ramipril_MFAT3 <- mean(matched32$high_fat,na.rm = TRUE)
ramipril_MFAT4 <- sd(matched32$high_fat,na.rm = TRUE)

# Calculate standardized mean difference (SMD) for age
sd_age_ramipril <- sqrt((ramipril_age2^2 + ramipril_age4^2) / 2)
smd_age_ramipril <- (ramipril_age1 - ramipril_age3) / sd_age_ramipril

# Calculate standardized mean difference (SMD) for BMI
sd_bmi_ramipril <- sqrt((ramipril_BMI2^2 + ramipril_BMI4^2) / 2)
smd_bmi_ramipril <- (ramipril_BMI1 - ramipril_BMI3) / sd_bmi_ramipril

# Calculate standardized mean difference (SMD) for ALT:
sd_ALT_ramipril <- sqrt((ramipril_ALT2^2 + ramipril_ALT4^2) / 2)
smd_ALT_ramipril <- (ramipril_ALT1 - ramipril_ALT3) / sd_ALT_ramipril

# Calculate standardized mean difference (SMD) for GLS:
sd_AST_ramipril <- sqrt((ramipril_AST2^2 + ramipril_AST4^2) / 2)
smd_AST_ramipril <- (ramipril_AST1 - ramipril_AST3) / sd_AST_ramipril

# Calculate standardized mean difference (SMD) for AST:
sd_GGT_ramipril <- sqrt((ramipril_GGT2^2 + ramipril_GGT4^2) / 2)
smd_GGT_ramipril <- (ramipril_GGT1 - ramipril_GGT3) / sd_GGT_ramipril

# Calculate standardized mean difference (SMD) for GGT:
sd_GLS_ramipril <- sqrt((ramipril_GLS2^2 + ramipril_GLS4^2) / 2)
smd_GLS_ramipril <- (ramipril_GLS1 - ramipril_GLS3) / sd_GLS_ramipril

# Calculate standardized mean difference (SMD) for age
sd_LFAT_ramipril <- sqrt((ramipril_LFAT2^2 + ramipril_LFAT4^2) / 2)
smd_LFAT_ramipril <- (ramipril_LFAT1 - ramipril_LFAT3) / sd_LFAT_ramipril

# Calculate standardized mean difference (SMD) for age
sd_MFAT_ramipril <- sqrt((ramipril_MFAT2^2 + ramipril_MFAT4^2) / 2)
smd_MFAT_ramipril <- (ramipril_MFAT1 - ramipril_MFAT3) / sd_MFAT_ramipril

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

ramipril_smd_diabetes1 <- prop_diabetes1 # Calculate standardized mean difference (SMD) for Diabetes being 1
ramipril_MASLD_yes <- sum(matched_ramipril$ramipril == 1 & matched_ramipril$MASLD == 1) # ramipril consuming MASLD patients : 
ramipril_MASLD_NO <- sum(matched_ramipril2$ramipril == 0 & matched_ramipril2$MASLD == 1)

# number of individuals consuming ramipril with MASLD
ramipril_individuals <- nrow(matched_ramipril)
ramipril_ind <- nrow(matched_ramipril2)

ramipril_MASLD_percent <- (ramipril_MASLD_yes / ramipril_individuals) * 100 # Calculate the percentage of individuals where both ramipril and MASLD:
non_ramipril_MASLD_percent <- (ramipril_MASLD_NO / ramipril_ind) * 100
ramipril_MASLD_sum <- sum(matched3$MASLD == 1, na.rm = TRUE) # Calculate the sum of binary values for Diabetes being 1
ramipril_MASLD_prop <- ramipril_MASLD_sum / length(matched_ramipril$MASLD) # Calculate the proportion of Diabetes being 1
ramipril_MASLD_n1 <- sum(!is.na(matched1$MASLD) & matched1$MASLD == 1) # Calculate the sample size for Diabetes being 1
ramipril_smd_MASLD <- ramipril_MASLD_prop # Calculate standardized mean difference (SMD) for Diabetes being 1
ramipril_MASH_yes <- sum(matched_ramipril$ramipril == 1 & matched_ramipril$MASH == 1) # ramipril consuming MASLD patients : 
ramipril_MASH_NO <- sum(matched_ramipril2$ramipril == 0 & matched_ramipril2$MASH == 1)
ramipril_MASH_percent <- (ramipril_MASH_yes / ramipril_individuals) * 100 # Calculate the percentage of individuals where both ramipril and MASH are equal to 1
non_ramipril_MASH_percent <- (ramipril_MASH_NO / ramipril_ind) * 100
ramipril_MASH_sum <- sum(matched3$MASH == 1, na.rm = TRUE) #  the sum of binary values for MASH being 1
ramipril_MASH_smd <- ramipril_MASH_sum / length(matched1$MASH) # Calculate the proportion of MASH being 1

# t test for ramipril

ramipril_t_test_MASLD <- t.test(matched3$ramipril , matched3$MASLD, na.rm = TRUE)
ramipril_t_test_MASH <- t.test(matched3$ramipril , matched3$MASH, na.rm = TRUE)
ramipril_t_test_Diabetes <- t.test(matched3$Diabetes, matched3$ramipril , na.rm = TRUE)
ramipril_t_test_Gender <- t.test(matched3$ramipril , matched3$Gender, na.rm = TRUE)
ramipril_t_test_ALT <- t.test(matched3$ramipril , matched3$ALT, na.rm = TRUE)
ramipril_t_test_ALT <- t.test(matched3$ramipril , matched3$GLS, na.rm = TRUE)

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

# ALT 
pioglitazone_ALT1 <- mean(matched21$ALT,na.rm = TRUE)
pioglitazone_ALT2 <- sd(matched21$ALT,na.rm = TRUE)
pioglitazone_ALT3 <- mean(matched22$ALT,na.rm = TRUE)
pioglitazone_ALT4 <- sd(matched22$ALT,na.rm = TRUE)

# AST:
pioglitazone_AST1 <- mean(matched21$AST,na.rm = TRUE)
pioglitazone_AST2 <- sd(matched21$AST,na.rm = TRUE)
pioglitazone_AST3 <- mean(matched22$AST,na.rm = TRUE)
pioglitazone_AST4 <- sd(matched22$AST,na.rm = TRUE)

# GGT 
pioglitazone_GGT1 <- mean(matched21$GGT,na.rm = TRUE)
pioglitazone_GGT2 <- sd(matched21$GGT,na.rm = TRUE)
pioglitazone_GGT3 <- mean(matched22$GGT,na.rm = TRUE)
pioglitazone_GGT4 <- sd(matched22$GGT,na.rm = TRUE)

# Glucose:
pioglitazone_GLS1 <- mean(matched21$GLS,na.rm = TRUE)
pioglitazone_GLS2 <- sd(matched21$GLS,na.rm = TRUE)
pioglitazone_GLS3 <- mean(matched22$GLS,na.rm = TRUE)
pioglitazone_GLS4 <- sd(matched22$GLS,na.rm = TRUE)

# Pooled standard deviation for age
pooled_sd_age <- sqrt((pioglitazone_age2^2 + pioglitazone_age4^2) / 2)
pioglitazone_smd_age <- (pioglitazone_age1 - pioglitazone_age3) / pooled_sd_age #Standardized mean difference for age

# Pooled standard deviation for BMI
pooled_sd_bmi <- sqrt((pioglitazone_BMI2^2 + pioglitazone_BMI4^2) / 2)
pioglitazone_smd_bmi <- (pioglitazone_BMI1 - pioglitazone_BMI3) / pooled_sd_bmi # Standardized mean difference for BMI

# Pooled standard deviation for ALT
pioglitazone_sd_ALT <- sqrt((pioglitazone_ALT2^2 + pioglitazone_ALT4^2) / 2)
pioglitazone_smd_ALT <- (pioglitazone_ALT1 - pioglitazone_ALT3) / pioglitazone_sd_ALT #Standardized mean difference for age

# Pooled standard deviation for Glucose
pioglitazone_sd_GLS <- sqrt((pioglitazone_GLS2^2 + pioglitazone_GLS4^2) / 2)
pioglitazone_smd_GLS <- (pioglitazone_GLS1 - pioglitazone_GLS3) / pioglitazone_sd_GLS #Standardized mean difference for age

# Pooled standard deviation for AST
pioglitazone_sd_AST <- sqrt((pioglitazone_AST2^2 + pioglitazone_AST4^2) / 2)
pioglitazone_smd_AST <- (pioglitazone_AST1 - pioglitazone_AST3) / pioglitazone_sd_AST #Standardized mean difference for age

# Pooled standard deviation for GGT
pioglitazone_sd_GGT <- sqrt((pioglitazone_GGT2^2 + pioglitazone_GGT4^2) / 2)
pioglitazone_smd_GGT <- (pioglitazone_GGT1 - pioglitazone_GGT3) / pioglitazone_sd_GGT #Standardized mean difference for age

# liver fat matched mean and SD for metformin and nonmetformin :

pioglitazone_mean_LFAT <- mean(matched21$low_fat)
pioglitazone_sd_LFAT <- sd(matched21$low_fat)
nonpioglitazone_mean_LFAT <- mean(matched22$low_fat)
nonpioglitazone_sd_LFAT <- sd(matched22$low_fat)

# Calculate pooled standard deviation for Low fat
pioglitazone_sd_LFAT <- sqrt((pioglitazone_sd_LFAT^2 + nonpioglitazone_sd_LFAT^2) / 2)
pioglitazone_smd_LFAT <- (pioglitazone_mean_LFAT - nonpioglitazone_mean_LFAT) / pioglitazone_sd_LFAT

# liver Low fat matched mean and SD for metformin and nonmetformin :

pioglitazone_mean_MFAT <- mean(matched21$high_fat)
pioglitazone_sd_MFAT <- sd(matched21$high_fat)
nonpioglitazone_mean_MFAT <- mean(matched22$high_fat)
nonpioglitazone_sd_MFAT <- sd(matched22$high_fat)

# Calculate pooled standard deviation for Low fat
pioglitazone_sd_MFAT <- sqrt((pioglitazone_sd_MFAT^2 + nonpioglitazone_sd_MFAT^2) / 2)
pioglitazone_smd_MFAT <- (pioglitazone_mean_MFAT - nonpioglitazone_mean_MFAT) / pioglitazone_sd_MFAT

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

pioglitazone_diabetes_yes <- sum(matched21$Diabetes == 1, na.rm = TRUE)
pioglitazone_diabetes_no <- sum(matched22$Diabetes == 0, na.rm = TRUE)
non_pioglitazone_diabetes_yes <- sum(matched22$Diabetes == 1, na.rm = TRUE)
non_pioglitazone_diabetes_no <- sum(matched22$Diabetes == 0, na.rm = TRUE)

# Calculate the sum of binary values for Diabetes being 1
pioglitazone_sum_diabetes1 <- sum(matched2$Diabetes, na.rm = TRUE)

# Calculate the proportion of Diabetes 
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
pioglitazone_MASH_yes <- sum(matched21$pioglitazone == 1 & matched21$MASH == 1)
pioglitazone_MASH_NO <- sum(matched22$pioglitazone == 0 & matched22$MASH == 1)

# Calculate the percentage of individuals where both pioglitazone and MASH are equal to 1
pioglitazone_MASH_percent <- (pioglitazone_MASH_yes / pioglitazone_individuals) * 100
non_pioglitazone_MASH_percent <- (pioglitazone_MASH_NO / pioglitazone_ind) * 100

#  the sum of binary values for MASH being 1
pioglitazone_MASH_sum <- sum(matched2$MASH == 1, na.rm = TRUE)
pioglitazone_MASH_smd <- pioglitazone_MASH_sum / length(matched2$MASH)

# t test for pioglitazone

pioglitazone_t_test_MASLD <- t.test(matched2$pioglitazone , matched2$MASLD, na.rm = TRUE)
pioglitazone_t_test_MASH <- t.test(matched2$pioglitazone , matched2$MASH, na.rm = TRUE)
pioglitazone_t_test_Diabetes <- t.test(matched2$pioglitazone , matched2$Diabetes, na.rm = TRUE)
pioglitazone_t_test_Gender <- t.test(matched2$pioglitazone , matched2$Gender, na.rm = TRUE)

# T-Test at Age baseline T-Test with/without metformin:
#calculating Mean and SD of Metformin with Age of baseline :

# Calculate mean and standard deviation of BMI
matched1_metformin <- subset(matched5, metformin == 1)
matched1_nonmetformin <- subset(matched5, metformin == 0)

# BMI matched mean and SD for metformin and nonmetformin :
metformin_mean_bmi <- mean(matched1_metformin$BMI)
metformin_sd_bmi <- sd(matched1_metformin$BMI)
nonmetformin_mean_bmi <- mean(matched1_nonmetformin$BMI)
nonmetformin_sd_bmi <- sd(matched1_nonmetformin$BMI)


# Calculate pooled standard deviation
metformin_sd_bmi <- sqrt((metformin_sd_bmi^2 + nonmetformin_sd_bmi^2) / 2)
smd_bmi <- (metformin_mean_bmi - nonmetformin_mean_bmi) / metformin_sd_bmi

# Glucose matched mean and SD for metformin and nonmetformin :
metformin_mean_GLS <- mean(matched1_metformin$GLS)
metformin_sd_GLS <- sd(matched1_metformin$GLS)
nonmetformin_mean_GLS <- mean(matched1_nonmetformin$GLS)
nonmetformin_sd_GLS <- sd(matched1_nonmetformin$GLS)

# Glucose standard deviation
metformin_sd_GLS <- sqrt((metformin_sd_GLS^2 + nonmetformin_sd_GLS^2) / 2)
smd_GLS <- (metformin_mean_GLS - nonmetformin_mean_GLS) / metformin_sd_GLS

# ALT matched mean and SD for metformin and nonmetformin :
metformin_mean_ALT <- mean(matched1_metformin$ALT)
metformin_sd_ALT <- sd(matched1_metformin$ALT)
nonmetformin_mean_ALT <- mean(matched1_nonmetformin$ALT)
nonmetformin_sd_ALT <- sd(matched1_nonmetformin$ALT)

# Calculate pooled standard deviation ALT:

metformin_sd_ALT <- sqrt((metformin_sd_ALT^2 + nonmetformin_sd_ALT^2) / 2)
smd_ALT <- (metformin_mean_ALT - nonmetformin_mean_ALT) / metformin_sd_ALT

# GGT matched mean and SD for metformin and nonmetformin :
metformin_mean_GGT <- mean(matched1_metformin$GGT)
metformin_sd_GGT <- sd(matched1_metformin$GGT)
nonmetformin_mean_GGT <- mean(matched1_nonmetformin$GGT)
nonmetformin_sd_GGT <- sd(matched1_nonmetformin$GGT)

# Calculate pooled standard deviation GGT :
metformin_sd_GGT <- sqrt((metformin_sd_GGT^2 + nonmetformin_sd_GGT^2) / 2)
smd_GGT <- (metformin_mean_GGT - nonmetformin_mean_GGT) / metformin_sd_GGT

# AST matched mean and SD for metformin and nonmetformin :

metformin_mean_AST <- mean(matched1_metformin$AST)
metformin_sd_AST <- sd(matched1_metformin$AST)
nonmetformin_mean_AST <- mean(matched1_nonmetformin$AST)
nonmetformin_sd_AST <- sd(matched1_nonmetformin$AST)

# Calculate pooled standard deviation
metformin_sd_AST <- sqrt((metformin_sd_AST^2 + nonmetformin_sd_AST^2) / 2)
smd_AST <- (metformin_mean_AST - nonmetformin_mean_AST) / metformin_sd_AST

# liver fat matched mean and SD for metformin and nonmetformin :

metformin_mean_LFAT <- mean(matched1_metformin$low_fat)
metformin_sd_LFAT <- sd(matched1_metformin$low_fat)
nonmetformin_mean_LFAT <- mean(matched1_nonmetformin$low_fat)
nonmetformin_sd_LFAT <- sd(matched1_nonmetformin$low_fat)

# Calculate pooled standard deviation for Low fat
metformin_sd_LFAT <- sqrt((metformin_sd_LFAT^2 + nonmetformin_sd_LFAT^2) / 2)
smd_LFAT <- (metformin_mean_LFAT - nonmetformin_mean_LFAT) / metformin_sd_LFAT

# liver Low fat matched mean and SD for metformin and nonmetformin :

metformin_mean_MFAT <- mean(matched1_metformin$high_fat)
metformin_sd_MFAT <- sd(matched1_metformin$high_fat)
nonmetformin_mean_MFAT <- mean(matched1_nonmetformin$high_fat)
nonmetformin_sd_MFAT <- sd(matched1_nonmetformin$high_fat)

# Calculate pooled standard deviation for Low fat
metformin_sd_MFAT <- sqrt((metformin_sd_MFAT^2 + nonmetformin_sd_MFAT^2) / 2)
smd_MFAT <- (metformin_mean_MFAT - nonmetformin_mean_MFAT) / metformin_sd_MFAT

# Age matched mean and SD :

metformin_mean_age <- mean(matched1_metformin$Age_AC)
matformin_sd_age <- sd(matched1_metformin$Age_AC)
nonmetformin_mean_age <- mean(matched1_nonmetformin$Age_AC)
nonmetformin_sd_age <- sd(matched1_nonmetformin$Age_AC)

# BMI :
nonmetformin_sd_age <- sqrt((matformin_sd_age^2 + nonmetformin_sd_age^2) / 2)
nonmetformin_smd_age <- (matformin_sd_age - nonmetformin_sd_age) / nonmetformin_sd_age

# Create a subset where metformin is 1 and MASH is 1
matched1_metformin_MASH <- subset(matched1_metformin, MASH == 1)
subset_count <- nrow(matched1_metformin_MASH)
total_count <- nrow(matched1_metformin)
percentage <- (subset_count / total_count) * 100

# Create a subset where metformin is 1 and MASH is 1
matched1_nonmetformin_MASH <- subset(matched1_nonmetformin, MASH == 1)
subset_count <- nrow(matched1_nonmetformin_MASH)
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
matched_smd_mafl <- (matched_prop_mafl1 - matched_prop_mafl2) / sqrt((matched1_n1_mafl_numeric * matched1_n2_mafl_numeric) / (matched1_n1_mafl_numeric + matched1_n2_mafl_numeric))

# Calculate proportions MASH
matched_prop_MASH1 <- sum(matched1$MASH == 0, na.rm = TRUE) / length(matched1$MASH)
matched_prop_MASH2 <- sum(matched1$MASH == 1, na.rm = TRUE) / length(matched1$MASH)

# Calculate sample sizes
matched1_n1_MASH <- sum(!is.na(matched1$MASH) & matched1$MASH == 0)
matched1_n2_MASH <- sum(!is.na(matched1$MASH) & matched1$MASH == 1)

# Convert sample sizes to numeric
matched1_n1_MASH_numeric <- as.numeric(matched1_n1_MASH)
matched1_n2_MASH_numeric <- as.numeric(matched1_n2_MASH)

# Calculate standardized mean difference (SMD) for MASH
matched_smd_MASH <- (matched_prop_MASH1 - matched_prop_MASH2) / sqrt((matched1_n1_MASH_numeric * matched1_n2_MASH_numeric) / (matched1_n1_MASH_numeric + matched1_n2_MASH_numeric))

# t test for metformin

matched1_t_test_MASLD <- t.test(matched1$metformin , matched1$MASLD, na.rm = TRUE)
matched1_t_test_MASH <- t.test(matched1$metformin , matched1$MASH, na.rm = TRUE)
matched1_t_test_Diabetes <- t.test(matched1$metformin , matched1$Diabetes, na.rm = TRUE)
matched1_t_test_Gender <- t.test(matched1$metformin , matched1$Gender, na.rm = TRUE)

# AlchoholicLiver dieases list: 

ramipril_alcoholicl_y <- sum(matched31$ramipril == 1 & matched31$alcoholicliver == 1)  
nonramipril_alcoholicl_N <- sum(matched32$ramipril == 0 & matched32$alcoholicliver == 1)

met_ind1 <- sum(!is.na(matched31$alcoholicliver) & matched31$alcoholicliver == 0)
met_ind2 <- sum(!is.na(matched3$alcoholicliver) & matched3$alcoholicliver == 1)

ramipril_alcoholicl_percent <- (ramipril_alcoholicl_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_alcoholicl_percent <- (nonramipril_alcoholicl_N / met_ind2) * 100

matched3_n1_ALCH <- sum(!is.na(matched31$alcoholicliver == 1)) # Calculate sample sizes
matched3_n2_ALCH <- sum(!is.na(matched32$alcoholicliver == 0))

# Convert sample sizes to numeric
matched3_n1_ALCH_numeric <- as.numeric(matched3_n1_ALCH)
matched3_n2_ALCH_numeric <- as.numeric(matched3_n2_ALCH)

# Calculate standardized mean difference (SMD) for toxicliver :
matched_smd_alch <- (ramipril_alcoholicl_y - nonramipril_alcoholicl_N) / sqrt((matched3_n1_ALCH_numeric * matched3_n2_ALCH_numeric) / (matched3_n1_ALCH_numeric + matched3_n2_ALCH_numeric))

# Toxic liver diseases list: 

ramipril_toxicliver_y <- sum(matched31$ramipril == 1 & matched31$toxicliver == 1) # ramipril consuming MASLD patients : 
nonramipril_toxicliver_N <- sum(matched32$ramipril == 0 & matched32$toxicliveriver == 1)

ramipril_toxicliver_percent <- (ramipril_toxicliver_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_toxicliver_percent <- (nonramipril_toxicliver_N / met_ind2) * 100

matched3_n1_tol <- sum(!is.na(matched3$toxicliver == 1)) # Calculate sample sizes
matched3_n2_tol <- sum(!is.na(matched3$toxicliver == 0))

# Convert sample sizes to numeric
matched3_n1_tol_numeric <- as.numeric(matched3_n1_tol)
matched3_n2_tol_numeric <- as.numeric(matched3_n2_tol)

# Calculate standardized mean difference (SMD) for Liverfailure :
matched_smd_tol <- (ramipril_toxicliver_y - nonramipril_toxicliver_N) / sqrt((matched3_n1_tol_numeric * matched3_n2_tol_numeric) / (matched3_n1_tol_numeric + matched3_n2_tol_numeric))

# liver failure diseases list: 

ramipril_liverfailure_y <- sum(matched31$ramipril == 1 & matched31$Liverfailure == 1) # ramipril consuming MASLD patients : 
nonramipril_liverfailure_N <- sum(matched32$ramipril == 0 & matched32$Liverfailure == 1)

ramipril_lvd_percent <- (ramipril_liverfailure_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_lvd_percent <- (nonramipril_liverfailure_N / met_ind2) * 100

matched3_n1_lvf <- sum(!is.na(matched3$Liverfailure == 1)) # Calculate sample sizes
matched3_n2_lvf <- sum(!is.na(matched3$Liverfailure == 0))

# Convert sample sizes to numeric
matched3_n1_lvf <- as.numeric(matched3_n1_lvf)
matched3_n2_lvf <- as.numeric(matched3_n2_lvf)

# Calculate standardized mean difference (SMD) for Liverfailure :
matched_smd_lvf <- (ramipril_liverfailure_y - nonramipril_liverfailure_N) / sqrt((matched3_n1_lvf * matched3_n2_lvf) / (matched3_n1_lvf + matched3_n2_lvf))

# ChronHepatitis diseases list: 

ramipril_chronc_y <- sum(matched31$ramipril == 1 & matched31$ChronHepatitis == 1) # ramipril consuming MASLD patients : 
nonramipril_chronc_N <- sum(matched32$ramipril == 0 & matched32$ChronHepatitis == 1)

ramipril_ChronHepatitis_percent <- (ramipril_chronc_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_ChronHepatitis_percent <- (nonramipril_chronc_N / met_ind2) * 100

matched3_n1_chronc <- sum(!is.na(matched31$ChronHepatitis == 0)) # Calculate sample sizes
matched3_n2_chronc <- sum(!is.na(matched32$ChronHepatitis == 1))

# Convert sample sizes to numeric
matched3_n1_cronc <- as.numeric(matched3_n1_chronc)
matched3_n2_cronc <- as.numeric(matched3_n2_chronc)

# Calculate standardized mean difference (SMD) for ChronHepatitis :
matched_smd_cronc <- (ramipril_chronc_y - nonramipril_chronc_N) / sqrt((matched3_n1_cronc * matched3_n2_cronc) / (matched3_n1_cronc + matched3_n2_cronc))

# Calculate proportions
ramipril_fbrcrh_y <- sum(matched31$ramipril == 1 & matched31$fibrosecirrho == 1)
nonramipril_fbrcrh_N <- sum(matched32$ramipril == 0 & matched32$fibrosecirrho == 1)

ramipril_fbrcrh_pr <- (ramipril_fbrcrh_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_fbrcrh_pr <- (nonramipril_fbrcrh_N / met_ind2) * 100

matched3_n1_fbrcrh <- sum(!is.na(matched31$fibrosecirrho == 0)) # Calculate sample sizes
matched3_n2_fbrcrh <- sum(!is.na(matched32$fibrosecirrho == 1))

# Calculate sample sizes
matched3_n1_fbrcrh <- as.numeric(matched3_n1_fbrcrh)
matched3_n2_fbrcrh <- as.numeric(matched3_n2_fbrcrh)

# Calculate standardized mean difference (SMD) for fibrosecirrhosis :
matched_smd_fbrcrh <- (ramipril_fbrcrh_y - nonramipril_fbrcrh_N) / sqrt((matched3_n1_fbrcrh * matched3_n2_fbrcrh) / (matched3_n1_fbrcrh + matched3_n2_fbrcrh))

# Calculate proportions
ramipril_inflam_y <- sum(matched31$ramipril == 1 & matched31$inflammliver == 1)
nonramipril_inflam_N <- sum(matched32$ramipril == 0 & matched32$inflammliver == 1)

ramipril_inflam_p <- (ramipril_inflam_y / met_ind1) * 100 # Calculate the percentage of individuals where both ramipril and MASLD are equal to 1
nonramipril_inflam_pr <- (nonramipril_inflam_N / met_ind2) * 100

# Calculate sample sizes
matched3_n1_inflamm <- sum(!is.na(matched31$inflammliver == 0)) 
matched3_n2_inflamm <- sum(!is.na(matched32$inflammliver == 1))

# Calculate sample sizes
matched3_n1_inflamm <- as.numeric(matched3_n1_inflamm)
matched3_n2_inflamm <- as.numeric(matched3_n2_inflamm)

# Calculate standardized mean difference (SMD) for fibrosecirrhosis :
matched_smd_inflamm <- (ramipril_inflam_y - nonramipril_inflam_N) / sqrt((matched3_n1_inflamm * matched3_n2_inflamm) / (matched3_n1_inflamm + matched3_n2_inflamm))

# Calculate proportions
ramipril_otherlD_y <- sum(matched31$ramipril == 1 & matched31$otherliverD == 1)
nonramipril_otherlD_N <- sum(matched32$ramipril == 0 & matched32$otherliverD == 1) 

# Calculate percentages
matched3_otherlD <- sum(ramipril_otherlD_y /met_ind1) * 100
matched3_otherlD2 <- sum(nonramipril_otherlD_N / met_ind2) * 100

# Calculate sample sizes
matched3_n1_otherlD <- sum(!is.na(matched3$otherliverD == 1)) 
matched3_n2_otherlD <- sum(!is.na(matched3$otherliverD == 0))

# Calculate sample sizes
matched3_n1_otherlD <- as.numeric(matched3_n1_otherlD)
matched3_n2_otherlD <- as.numeric(matched3_n2_otherlD)

# Calculate standardized mean difference (SMD) for fibrosecirrhosis :
matched_smd_otherlD <- (ramipril_otherlD_y - nonramipril_otherlD_N) / sqrt((matched3_n1_otherlD * matched3_n2_otherlD) / (matched3_n1_otherlD + matched3_n2_otherlD))

# Calculate proportions
ramipril_Unclassified_y <- sum(matched31$ramipril == 1 & matched31$UnclassLiverd == 1)
nonramipril_Unclassified_N <- sum(matched32$ramipril == 0 & matched32$UnclassLiverd == 1) 

# Calculate percentages
matched3_Unclassified <- sum(matched31$UnclassLiverd / met_ind1) * 100
matched3_Unclassified2 <- sum(matched32$UnclassLiverd / met_ind2) * 100

# Calculate sample sizes
matched3_n1_Unclassified <- sum(!is.na(matched3$UnclassLiverd == 0)) 
matched3_n2_Unclassified <- sum(!is.na(matched3$UnclassLiverd == 1))

# Calculate sample sizes
matched3_n1_Unclassified <- as.numeric(matched3_n1_Unclassified)
matched3_n2_Unclassified <- as.numeric(matched3_n2_Unclassified)

# Calculate standardized mean difference (SMD) for fibrosecirrhosis :
matched_smd_Unclassified <- (ramipril_Unclassified_y - nonramipril_Unclassified_N) / sqrt((matched3_n1_Unclassified * matched3_n2_Unclassified) / (matched3_n1_Unclassified + matched3_n2_Unclassified))
