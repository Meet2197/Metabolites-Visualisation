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

# Identify rows that are not related to metformin and create negative control dataframe
non_metformin <- medication[!medication$eid %in% metformin_table$eid, ]

# Create data table for negative control
non_metformin = non_metformin[!duplicated(non_metformin), ]
non_metformin <- data.table(non_metformin)


#MASLD icd code K760 to create new df. 

MASLD <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K760'))
MASLD <- as.data.frame(MASLD[!duplicated(MASLD$eid_1), ])
MASLD$Diagnosis <- 'K760'
MASLD <- select(MASLD, c(eid_1, Diagnosis))
setnames(MASLD,"eid_1","eid")

# NON MASLD dataframe : 

NON_MASLD <- subset(hesin_diag, !startsWith(as.character(diag_icd10), 'K760'))
NON_MASLD <- as.data.frame(NON_MASLD[!duplicated(NON_MASLD$eid_1), ])

# NASH dataframe K758 ICD CODE: 

NASH <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K758'))
NASH <- as.data.frame(NASH[!duplicated(NASH$eid_1), ])
NASH$Diagnosis <- 'K758'
NASH <- select(NASH, c(eid_1, Diagnosis))
setnames(NASH,"eid_1","eid")

# NON_NASH :

NON_NASH <- subset(hesin_diag, !startsWith(as.character(diag_icd10), 'K758'))
NON_NASH <- as.data.frame(NON_NASH[!duplicated(NON_NASH$eid_1), ])

# MASLD with metformin: 

common7 <- intersect(metformin_df$eid, MASLD$eid_1)
metformin_common <- metformin_df[metformin_df$eid %in% common7, ]
MASLD_common <- MASLD[MASLD$eid_1 %in% common7, ]
metformin_MASLD <- merge(metformin_common, MASLD_common, 
                        by.x = "eid", by.y = "eid_1", suffixes = c('_metformin', '_MASLD'))
metformin_MASLD <- metformin_MASLD %>% 
  select(eid, Diagnosis)

# MASLD metformin percentage calculation against metformin df:

MASLD_metformin_percent <- length(intersect(metformin_MASLD$eid, metformin_df$eid)) / length(metformin_df$eid) * 100


# T test analysis with MASLD and Metformin 
metformin_table$Metformin<-1
MASLD$MASLD<-1

Masld_t_test_metformin<-merge(covariates, metformin_table, by.x="eid", all.x=TRUE)
Masld_t_test_metformin<-merge(Masld_t_test_metformin, MASLD, by.x="eid", by.y ="eid")
Masld_t_test_metformin$MASLD<-replace_na(Masld_t_test_metformin$MASLD, 0)
Masld_t_test_metformin$Metformin<-replace_na(Masld_t_test_metformin$Metformin, 0)
MASLD_t_test <- t.test(Masld_t_test_metformin$MASLD, Masld_t_test_metformin$Metformin)

Masld_t_test_metformin$MASLD <- as.numeric(Masld_t_test_metformin$MASLD) - 1
Masld_t_test_metformin$Metformin <- as.numeric(Masld_t_test_metformin$Metformin) - 1
MASLD$eid_1 <- as.numeric(MASLD$eid_1) - 1
metformin_table$eid <- as.numeric(metformin_table$eid) - 1
setnames(MASLD,"eid_1","eid")

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

# T-Test for MASLD:

metformin_MASLD$eid <- as.numeric(metformin_MASLD$eid) - 1
MASLD_non_metformin$eid <- as.numeric(MASLD_non_metformin$eid) - 1

# T test for Non metformin :

MASLD$MASLD<-1
non_metformin$nonmetformin<-1

Masld_t_test_nonmetformin<-merge(covariates, non_metformin, by.x="eid", all.x=TRUE)
Masld_t_test_nonmetformin<-merge(Masld_t_test_nonmetformin, MASLD, by.x="eid", all.x=TRUE)
Masld_t_test_nonmetformin$MASLD<-replace_na(Masld_t_test_nonmetformin$MASLD, 0)
Masld_t_test_nonmetformin$nonmetformin<-replace_na(Masld_t_test_nonmetformin$nonmetformin, 0)
MASLD_t_test_nonmetformin <- t.test(Masld_t_test_nonmetformin$MASLD, Masld_t_test_nonmetformin$nonmetformin)

Masld_t_test_nonmetformin$MASLD <- as.numeric(Masld_t_test_nonmetformin$MASLD) - 1
Masld_t_test_nonmetformin$nonmetformin <- as.numeric(Masld_t_test_nonmetformin$nonmetformin) - 1

# Mean and standard deviation with/without metformin dataframes:

mean_MASLD_metformin <- mean(metformin_MASLD$eid, na.rm = TRUE)
mean_MASLD_non_metformin <- mean(MASLD_non_metformin$eid, na.rm = TRUE)
sd_MASLD_metformin<- sd(metformin_MASLD$eid, na.rm = TRUE)
sd_MASLD_non_metformin <- sd(MASLD_non_metformin$eid, na.rm = TRUE)

# chanigng data with numeric format:  

mean_MASLD_metformin <- as.numeric(mean_MASLD_metformin)
mean_MASLD_non_metformin <- as.numeric(mean_MASLD_non_metformin)
sd_MASLD_metformin <- as.numeric(sd_MASLD_metformin)
sd_MASLD_non_metformin <- as.numeric(sd_MASLD_non_metformin)

# standardized mean difference for MASLD ICD 10 with/without metformin:

smd_MASLD <- (mean_MASLD_metformin - mean_MASLD_non_metformin) / sqrt((sd_MASLD_metformin^2 + sd_MASLD_non_metformin^2) / 2)

# NASH T-Test with metformin for P value generation:

NASH$NASH <- 1
metformin$metformin <- 1

NASH_t_metformin<-merge(covariates, metformin, by.x="eid", all.x=TRUE)
NASH_t_metformin<-merge(NASH_t_metformin, NASH, by.x="eid", all.x=TRUE)
NASH_t_metformin$NASH<-replace_na(NASH_t_metformin$NASH, 0)
NASH_t_metformin$metformin<-replace_na(NASH_t_metformin$metformin, 0)
NASH_t_metformin$NASH <- as.numeric(NASH_t_metformin$NASH)
NASH_t_metformin$metformin <- as.numeric(NASH_t_metformin$metformin)

NASH_t_metformin <- t.test(NASH_t_metformin$NASH, NASH_t_metformin$metformin)

setnames(NASH,"eid_1","eid")

# NASH with metformin: 

common9 <- intersect(metformin_df$eid, NASH$eid_1)
metformin_common <- metformin_df[metformin_df$eid %in% common9, ]
NASH_common <- NASH[NASH$eid_1 %in% common9, ]
NASH_metformin <- merge(metformin_common, NASH_common, 
                        by.x = "eid", by.y = "eid_1", suffixes = c('_metformin', '_NASH'))
NASH_metformin <- NASH_metformin %>% 
  select(eid, Diagnosis)

# NASH metformin percentage calculation: 

NASH_metformin_percent <- length(intersect(NASH_metformin$eid, metformin_df$eid)) / length(metformin_df$eid) * 100

# NASH with Non metformin : 

common10 <- intersect(non_metformin$eid, NASH$eid_1)
non_metformin_common <- non_metformin[non_metformin$eid %in% common10, ]
NASH_common <- NASH[NASH$eid_1 %in% common10, ]
non_metformin_common$eid <- as.integer(non_metformin_common$eid)
NASH_non_metformin <- merge(non_metformin_common, NASH_common, 
                                by.x = "eid", by.y = "eid_1", suffixes = c('_Non_metformin', '_NASLH'))
NASH_non_metformin <- NASH_non_metformin %>% 
  select(eid, Diagnosis)

# NASH with metformin percentage calculation : 

NASH_nonmetformin_percent <- length(intersect(NASH_non_metformin$eid, non_metformin$eid)) / length(non_metformin$eid) * 100

# NASH numeric format changes : 

NASH_metformin$eid <- as.numeric(NASH_metformin$eid) - 1
NASH_non_metformin$eid <- as.numeric(NASH_non_metformin$eid) - 1

# NASH T-Test with non meyformin for P value generation:

NASH$NASH<-1
MASLD$MASLD<-1
non_metformin$nonmetformin<-1
metformin$metformin<-1

merged_t_df<-merge(covariates, non_metformin, by.x="eid", all.x=TRUE)
merged_t_df<-merge(merged_t_df, MASLD, by.x="eid", by.y = "eid") %>%
  select('eid','BMI','Age_at_Obs.','Gender','nonmetformin','MASLD')

merged_t_df<-merge(merged_t_df, NASH, by.x="eid", by.y="eid") %>%
  select('eid','BMI','Age_at_Obs.','Gender','nonmetformin','MASLD','NASH')
merged_t_df<-merge(merged_t_df, metformin, by.x="eid", by.y="eid") %>%
  select('eid','BMI','Age_at_Obs.','Gender','nonmetformin','MASLD','NASH','metformin')

merged_t_df$MASLD<-replace_na(merged_t_df$MASLD, 0)
merged_t_df$metformin<-replace_na(merged_t_df$metformin, 0)
merged_t_df$nonmetformin<-replace_na(merged_t_df$nonmetformin, 0)
merged_t_df$NASH<-replace_na(merged_t_df$NASH, 0)

Age_MASLD_t_metformin <- t.test(merged_t_df$MASLD, merged_t_df$metformin)

merged_t_df$MASLD <- as.numeric(merged_t_df$MASLD)
merged_t_df$nonmetformin <- as.numeric(merged_t_df$nonmetformin)


mean_NASH_metformin <- mean(NASH_metformin$eid, na.rm = TRUE)
mean_NASH_non_metformin <- mean(NASH_non_metformin$eid, na.rm = TRUE)

# Numeric conversion for NASH with Metformin and non metformin: 

mean_NASH_metformin <- as.numeric(mean_NASH_metformin)
mean_NASH_non_metformin <- as.numeric(mean_NASH_non_metformin)
sd_NASH_metformin <- as.numeric(sd_NASH_metformin)
sd_NASH_non_metformin <- as.numeric(sd_NASH_non_metformin)


sd_NASH_metformin<- sd(NASH_metformin$eid, na.rm = TRUE)
sd_NASH_non_metformin <- sd(NASH_non_metformin$eid, na.rm = TRUE)


# standardized mean difference for NASH from ICD 10 with/without metformin:

smd_NASH <- (mean_NASH_metformin - mean_NASH_non_metformin) / sqrt((sd_NASH_metformin^2 + sd_NASH_non_metformin^2) / 2)

# Metformin with Age of baseline characteristics:

common11 <- intersect(metformin_df$eid, baseline_df$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common11, ]
baseline_df_common  <- baseline_df[baseline_df$eid %in% common11, ]
age_metformin <- merge(metformin_df_common, baseline_df_common, 
                                by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_Age'))
age_metformin <- age_metformin %>% 
  select(eid, Age_AC)

# T-Test at Age baseline T-Test with/without metformin:


NASH_t_nonmetformin<-merge(covariates, non_metformin, by.x="eid", all.x=TRUE)
NASH_t_nonmetformin<-merge(NASH_t_nonmetformin, NASH, by.x="eid", all.x=TRUE)
NASH_t_nonmetformin$NASH<-replace_na(NASH_t_nonmetformin$NASH, 0)
NASH_t_nonmetformin$non_metformin<-replace_na(NASH_t_nonmetformin$non_metformin, 0)
NASH_t_nonmetformin <- t.test(NASH_t_nonmetformin$NASH, NASH_t_nonmetformin$non_metformin)

NASH_t_nonmetformin$NASH <- as.numeric(NASH_t_nonmetformin$NASH)
NASH_t_nonmetformin$non_metformin <- as.numeric(NASH_t_nonmetformin$nonmetformin)


#calculating Mean and SD of Metformin with Age of baseline :

mean_metformin <- mean(age_metformin$Age_AC)
sd_metformin <- sd(age_metformin$Age_AC)

# non-Metformin with Age at baseline :

common12 <- intersect(non_metformin$eid, baseline_df$eid)
non_metformin_common <- non_metformin[non_metformin$eid %in% common12, ]
baseline_Age_common  <- baseline_df[baseline_df$eid %in% common12, ]
non_metformin_common$eid <- as.integer(non_metformin_common$eid)
baseline_non_metformin <- merge(non_metformin_common, baseline_Age_common,
                                by.x = "eid", by.y = "eid", suffixes = c('_non_metformin', '_Age'))
AgeAc_nonmetformin <- baseline_non_metformin %>% 
  select(eid, Age_AC)

# Mean and SD calculation for non metformin with Age at baseline :  

mean_non_metformin <- mean(AgeAc_nonmetformin$Age_AC)
sd_non_metformin <- sd(AgeAc_nonmetformin$Age_AC)

# Numeric format changes for Age at baseline with/without Metformin: 

age_metformin$eid <- as.numeric(age_metformin$eid) - 1
AgeAc_nonmetformin$eid <- as.numeric(AgeAc_nonmetformin$eid) - 1

#  standardized mean difference for Age at baseline with/without metformin:

age_smd <- (mean_metformin - mean_non_metformin) / sqrt((sd_metformin^2 + sd_non_metformin^2) / 2)

# Metformin with BMI :

common13 <- intersect(metformin_df$eid, baseline_df$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common13, ]
baseline_df_common  <- baseline_df[baseline_df$eid %in% common13, ]
baseline_metformin <- merge(metformin_df_common, baseline_df_common, 
                            by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_BMI'))
BMI_metformin <- baseline_metformin %>% 
  select(eid, BMI)

# Mean and SD calculation for non metformin with Age at baseline :  

mean_metformin <- mean(BMI_metformin$BMI, na.rm = TRUE)
sd_metformin <- sd(BMI_metformin$BMI, na.rm = TRUE)

# non-Metformin with BMI :

common14 <- intersect(non_metformin$eid, baseline_df$eid)
non_metformin_common <- non_metformin[non_metformin$eid %in% common14, ]
baseline_Age_common  <- baseline_df[baseline_df$eid %in% common14, ]
non_metformin_common$eid <- as.integer(non_metformin_common$eid)
baseline_non_metformin <- merge(non_metformin_common, baseline_Age_common,
                                by.x = "eid", by.y = "eid", suffixes = c('_non_metformin', '_BMI'))
BMI_nonmetformin <- baseline_non_metformin %>% 
  select(eid, BMI)

# Mean and SD calculation at BMI baseline non metformin :

mean_non_metformin <- mean(BMI_nonmetformin$BMI, na.rm = TRUE)
sd_nonmetformin <- sd(BMI_nonmetformin$BMI, na.rm = TRUE)

# Calculate standardized mean difference :
smd_before_psm <- (mean_metformin - mean_non_metformin) / sqrt((sd_metformin^2 + sd_nonmetformin^2) / 2)

BMI_metformin$eid <- as.numeric(BMI_metformin$eid) - 1
BMI_nonmetformin$eid <- as.numeric(BMI_nonmetformin$eid) - 1

# T-Test at BMI baseline T-Test with/without metformin::

BMI_t_test <- t.test(BMI_metformin$eid, BMI_nonmetformin$eid)


# Gender based Calculation of Metformin(Female)

baseline_df$Gender <- str_replace_all(baseline_df$Gender, c("0" = "Female", "1" = "Male"))
baseline_female <- baseline_df[baseline_df$Gender == "Female", ]
common15 <- intersect(metformin_df$eid, baseline_female$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common15, ]
baseline_female_common  <- baseline_female[baseline_female$eid %in% common15, ]
female_metformin <- merge(metformin_df_common, baseline_female_common, 
                            by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_Female'))
female_metformin <- female_metformin %>% 
  select(eid, Gender)

# Female gender at Baseline with metformin percentage calculation:

female_metformin_percent <- length(intersect(female_metformin$eid, metformin_df$eid)) / length(metformin_df$eid) * 100

# Gender based Calculation of NonMetformin(Female)

baseline_female <- baseline_df[baseline_df$Gender == "Female", ]
common16 <- intersect(non_metformin$eid, baseline_female$eid)
nonmetformin_df_common <- non_metformin[non_metformin$eid %in% common16, ]
baseline_female_common  <- baseline_female[baseline_female$eid %in% common16, ]
nonmetformin_df_common$eid <- as.integer(nonmetformin_df_common$eid)
female_nonmetformin <- merge(nonmetformin_df_common, baseline_female_common, 
                          by.x = "eid", by.y = "eid", suffixes = c('_nonmetformin', '_Female'))
female_nonmetformin <- female_nonmetformin %>% 
  select(eid, Gender)

# Female gender at Baseline with Non-metformin percentage calculation:

female_nonmetformin_percent <- length(intersect(female_nonmetformin$eid, non_metformin$eid)) / length(non_metformin$eid) * 100

female_nonmetformin$eid <- as.numeric(female_nonmetformin$eid) - 1
female_metformin$eid <- as.numeric(female_metformin$eid) - 1

# T-Test at Female gender baseline T-Test with/without metformin:

female_gender_t_test <- t.test(female_metformin$eid, female_nonmetformin$eid)

mean_female_metformin <- mean(female_metformin$eid, na.rm = TRUE)
mean_female_nonmetformin <- mean(female_nonmetformin$eid, na.rm = TRUE)
sd_female_metformin<- sd(female_metformin$eid, na.rm = TRUE)
sd_female_non_metformin <- sd(female_nonmetformin$eid, na.rm = TRUE)

# Numeric format changes for Female Gender with/without Metformin

mean_female_metformin <- as.numeric(mean_female_metformin)
mean_female_nonmetformin <- as.numeric(mean_female_nonmetformin)
sd_female_metformin <- as.numeric(sd_female_metformin)
sd_female_non_metformin <- as.numeric(sd_female_non_metformin)

# standardized mean difference for Female gender at baseline with/without metformin:

smd_female <- (mean_female_metformin - mean_female_nonmetformin) / sqrt((sd_female_metformin^2 + sd_female_non_metformin^2) / 2)

# Gender based Calculation of Metformin(Male)

baseline_df$Gender <- str_replace_all(baseline_df$Gender, c("0" = "Female", "1" = "Male"))
baseline_male <- baseline_df[baseline_df$Gender == "Male", ]
common17 <- intersect(metformin_df$eid, baseline_male$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common17, ]
baseline_male_common  <- baseline_male[baseline_male$eid %in% common17, ]
male_metformin <- merge(metformin_df_common, baseline_male_common, 
                          by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_Male'))
male_metformin <- male_metformin %>% 
  select(eid, Gender)

# male gender at Baseline with metformin percentage calculation:

male_metformin_percent <- length(intersect(male_metformin$eid, metformin_df$eid)) / length(metformin_df$eid) * 100

# Gender based Calculation of NonMetformin(Male)

common18 <- intersect(non_metformin$eid, baseline_male$eid)
nonmetformin_df_common <- non_metformin[non_metformin$eid %in% common18, ]
baseline_male_common  <- baseline_male[baseline_male$eid %in% common18, ]
nonmetformin_df_common$eid <- as.integer(nonmetformin_df_common$eid)
male_nonmetformin <- merge(nonmetformin_df_common, baseline_male_common, 
                             by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_Male'))
male_nonmetformin <- male_nonmetformin %>% 
  select(eid, Gender)
male_nonmetformin_percent <- length(intersect(male_nonmetformin$eid, baseline_df$eid)) / length(baseline_df$eid) * 100

male_nonmetformin$eid <- as.numeric(male_nonmetformin$eid) - 1
male_metformin$eid <- as.numeric(male_metformin$eid) - 1

# T-Test at Male Gender baseline T-Test with/without metformin:

male_gender_t_test <- t.test(male_nonmetformin$eid, male_metformin$eid)

mean_male_metformin <- mean(male_metformin$eid, na.rm = TRUE)
mean_male_nonmetformin <- mean(male_nonmetformin$eid, na.rm = TRUE)
sd_male_metformin<- sd(male_metformin$eid, na.rm = TRUE)
sd_male_non_metformin <- sd(male_nonmetformin$eid, na.rm = TRUE)

# Numeric format changes for Male Gender with/without Metformin

mean_male_metformin <- as.numeric(mean_male_metformin)
mean_male_nonmetformin <- as.numeric(mean_male_nonmetformin)
sd_male_metformin <- as.numeric(sd_male_metformin)
sd_male_non_metformin <- as.numeric(sd_male_non_metformin)

# standardized mean difference for Male Gender at baseline with/without metformin:

smd_male <- (mean_male_metformin - mean_male_nonmetformin) / sqrt((sd_male_metformin^2 + sd_male_non_metformin^2) / 2)

# Diabetes Calculation of Metformin

baseline_df$Diabetes <- str_replace_all(baseline_df$Diabetes, c("0" = "No", "1" = "Yes"))
diabetes_no <- baseline_df[baseline_df$Diabetes == "No", ]
common19 <- intersect(metformin_df$eid, diabetes_no$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common19, ]
diabetes_no_common  <- diabetes_no[diabetes_no$eid %in% common19, ]
diabetes_no_metformin <- merge(metformin_df_common, diabetes_no_common, 
                          by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_no'))
diabetes_no_metformin <- diabetes_no_metformin %>% 
  select(eid, Diabetes)

# Diabetes No at Baseline with metformin percentage calculation:

no_metformin_diabetes <- length(intersect(diabetes_no_metformin$eid, metformin_df$eid)) / length(metformin_df$eid) * 100

# Diabetes based Calculation of NonMetformin

common20 <- intersect(non_metformin$eid, diabetes_no$eid)
nonmetformin_df_common <- non_metformin[non_metformin$eid %in% common20, ]
diabetes_no_common  <- diabetes_no[diabetes_no$eid %in% common20, ]
nonmetformin_df_common$eid <- as.integer(nonmetformin_df_common$eid)
diabetes_no_nonmetformin <- merge(nonmetformin_df_common, diabetes_no_common, 
                             by.x = "eid", by.y = "eid", suffixes = c('_nonmetformin', '_no'))
diabetes_no_nonmetformin <- diabetes_no_nonmetformin %>% 
  select(eid, Diabetes)

# Diabetes No at Baseline with non-metformin percentage calculation:

no_nonmetformin_diabetes <- length(intersect(diabetes_no_nonmetformin$eid, non_metformin$eid)) / length(non_metformin$eid) * 100

diabetes_no_metformin$eid <- as.numeric(diabetes_no_metformin$eid) - 1
diabetes_no_nonmetformin$eid <- as.numeric(diabetes_no_nonmetformin$eid) - 1

# T-Test at No Diabetes of baseline T-Test with/without metformin:

no_diabetes_t_test <- t.test(diabetes_no_metformin$eid, diabetes_no_nonmetformin$eid)

mean_no_metformin <- mean(diabetes_no_metformin$eid, na.rm = TRUE)
mean_no_nonmetformin <- mean(diabetes_no_nonmetformin$eid, na.rm = TRUE)
sd_no_metformin<- sd(diabetes_no_metformin$eid, na.rm = TRUE)
sd_no_non_metformin <- sd(diabetes_no_nonmetformin$eid, na.rm = TRUE)

# Numeric format changes for No Diabetes with/without Metformin

mean_no_metformin <- as.numeric(mean_no_metformin)
mean_no_nonmetformin <- as.numeric(mean_no_nonmetformin)
sd_no_metformin <- as.numeric(sd_no_metformin)
sd_no_non_metformin <- as.numeric(sd_no_non_metformin)

# standardized mean difference for No of Diabetes at baseline with/without metformin:

smd_no <- (mean_no_metformin - mean_no_nonmetformin) / sqrt((sd_no_metformin^2 + sd_no_non_metformin^2) / 2)


# Diabetes based Calculation of Metformin

baseline_yes <- baseline_df[baseline_df$Diabetes == "Yes", ]
common21 <- intersect(metformin_df$eid, baseline_yes$eid)
metformin_df_common <- metformin_df[metformin_df$eid %in% common21, ]
baseline_yes_common  <- baseline_yes[baseline_yes$eid %in% common21, ]
yes_metformin <- merge(metformin_df_common, baseline_yes_common, 
                        by.x = "eid", by.y = "eid", suffixes = c('_metformin', '_Yes'))
yes_metformin <- yes_metformin %>% 
  select(eid, Diabetes)

# Diabetes Yes at Baseline with metformin percentage calculation:

yes_metformin_diabetes <- length(intersect(yes_metformin$eid, metformin_df$eid)) / length(metformin_df$eid) * 100


# Diabetes based Calculation of NonMetformin
baseline_yes <- baseline_df[baseline_df$Diabetes == "Yes", ]
common22 <- intersect(non_metformin$eid, baseline_yes$eid)
nonmetformin_df_common <- non_metformin[non_metformin$eid %in% common22, ]
baseline_yes_common  <- baseline_yes[baseline_yes$eid %in% common22, ]
nonmetformin_df_common$eid <- as.integer(nonmetformin_df_common$eid)
yes_nonmetformin <- merge(nonmetformin_df_common, baseline_yes_common, 
                           by.x = "eid", by.y = "eid", suffixes = c('_non_metformin', '_Yes'))
yes_nonmetformin <- yes_nonmetformin %>% 
  select(eid, Diabetes)

# Diabetes Yes at Baseline with non-metformin percentage calculation:

yes_nonmetformin_diabetes <- length(intersect(yes_nonmetformin$eid, non_metformin$eid)) / length(non_metformin$eid) * 100

yes_nonmetformin$eid <- as.numeric(yes_nonmetformin$eid) - 1
yes_metformin$eid <- as.numeric(yes_metformin$eid) - 1

# T-Test for Diabetes Yes : 

yes_metformin_t_test_result <- t.test(yes_nonmetformin$eid, yes_metformin$eid)

# Mean and SD calculation for non metformin with Yes at Diabetes baseline :  

mean_yes_metformin <- mean(yes_metformin$eid, na.rm = TRUE)
mean_yes_nonmetformin <- mean(yes_nonmetformin$eid, na.rm = TRUE)
sd_yes_metformin<- sd(yes_metformin$eid, na.rm = TRUE)
sd_yes_non_metformin <- sd(yes_nonmetformin$eid, na.rm = TRUE)

# Numeric format changes for Yes Diabetes with/without Metformin

mean_yes_metformin <- as.numeric(mean_yes_metformin)
mean_yes_nonmetformin <- as.numeric(mean_yes_nonmetformin)
sd_yes_metformin <- as.numeric(sd_yes_metformin)
sd_yes_non_metformin <- as.numeric(sd_yes_non_metformin)

# standardized mean difference for yes of Diabetes at baseline with/without metformin: 

smd_yes <- (mean_yes_metformin - mean_yes_nonmetformin) / sqrt((sd_yes_metformin^2 + sd_yes_non_metformin^2) / 2)
