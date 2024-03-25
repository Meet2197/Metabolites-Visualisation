# Death data with baseline examination: # Read the required data files :

requirement <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/baselineRequire.csv') %>% select('eid', 'X6138.0.0', 'X2443.0.0')
covariates2 <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates2.csv') %>% select('eid', 'X21003.0.0', 'X20116.0.0', 'X20117.0.0')
metabolites_baseline <- metabolites %>% select('Eid', 'BMI', 'Sex', 'Ethnicity')

# Read baseline data file :

baseline_df <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/baseline.csv')
hypertension <- read.csv(file = 'C:/Users/User/Desktop/Data/hypertension_baseline.csv')
nutritrion <- read.csv(file = 'C:/Users/User/Desktop/Data/nutritrion.csv')
init_exam <- read.csv(file = 'C:/Users/User/Desktop/Data/init.csv')

# Rename column names :

baseline_df <- as.data.frame(baseline_df[!duplicated(baseline_df$eid),])

# Read another data file :

# baseline <- fread("C:/Users/User/Desktop/Data/ukb52200.csv", select=c("eid","21003-0.0","31-0.0","21001-0.0", "20116-0.0", "20117-0.0", "41270-0.0", "6138-0.0", "2443-0.0"))
# hypertension_baseline_df <- fread("C:/Users/User/Desktop/Data/ukb52200.csv" ,select=c("eid", "21003-0.0", "31-0.0","20116-0.0", "20117-0.0", "41270-0.0", "6138-0.0", "2443-0.0", "4080-0.0","4079-0.0"))
# nutritrion <- fread("C:/Users/User/Desktop/Data/ukb52200.csv", select=c("eid","20089-0.0","1558-0.0","100580-0.0","1548-0.0","1309-0.0","1369-0.0","6164-0.0","1349-0.0","6144-0.0","1160-0.0","20090-0.0","4537-0.0","1930-0.0"))
# covariates <- fread("C:/Users/User/Desktop/Data/ukb52200.csv", select=c("eid","30620-0.0","30650-0.0","30730-0.0"))
# init_exam <- fread("C:/Users/User/Desktop/Data/ukb52200.csv", select=c("eid","53-0.0"))

# setnames for hypertension baseline_df:

setnames(baseline_df, old = colnames(baseline_df), new = c('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking','Diagnosis', 'Qualifications', 'Diabetes'))
setnames(hypertension, old = colnames(hypertension), new = c('eid','Age_AC', 'Smoking', 'Drinking','Diagnosis', 'Qualifications', 'Diabetes', 'Systolic_BP','Diastolic_BP'))
Covariates2<- fread("ukb52200.csv", select=c("eid", "21000-0.0","21001-0.0"))
setnames(nutritrion, old = colnames(nutritrion), new = c('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking','Diagnosis', 'Qualifications', 'Diabetes', 'meals', 'spz_diet', 'Alcohol_int_10y', 'Alcohol_int_fr','Alcohol_cons'))
setnames(init_exam, old = colnames(init_exam), new = c('eid', 'init_exam'))

# Create a new column for Systolic_BP categories # icd i11-15 observation 

hypertension_fr <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'I10') | startsWith(as.character(diag_icd10), 'I11') | startsWith(as.character(diag_icd10), 'I12') | startsWith(as.character(diag_icd10), 'I13')
                          | startsWith(as.character(diag_icd10), 'I14') | startsWith(as.character(diag_icd10), 'I15')) # ICD code filter
hypertension_fr = hypertension_fr[!duplicated(hypertension_fr$eid), ]

# K 70/77 merge with Liver disease df.

hesin_diag$diag_icd10<-replace_na(hesin_diag$diag_icd10, 0)
patterns <- c("K70", "K71", "K72", "K73", "K74", "K75", "K76", "K77")
regex_pattern <- paste0("^", paste(patterns, collapse = "|")) # regular expression pattern to match any of the patterns defined
liverdisease <- hesin_diag %>% filter(substr(diag_icd10, 1, 3) >= 'K70' & substr(diag_icd10, 1, 3) <= 'K77')

# Filter the dataframe based on the regex
#medication search 

hypertension$Systolic_Category <- cut(hypertension$Systolic_BP, breaks = c(-Inf, 120, 129, 139, 180, Inf), labels = c("Normal", "Elevated", "Stage 1", "Stage 2", "Stage 3"), include.lowest = TRUE)

# Create a new column for Diastolic_BP categories

hypertension$Diastolic_Category <- cut(hypertension$Diastolic_BP, breaks = c(-Inf, 80, 89, 120, Inf), labels = c("Normal", "Stage 1", "Stage 2", "Stage 3"), include.lowest = TRUE)

# Create the Hypertension column based on the conditions
hypertension$Hypertension <- ifelse(hypertension$Systolic_BP >= 140 | hypertension$Diastolic_BP >= 90, 1, 0)

# Rename column names :

setnames(covariates2, old = colnames(covariates2), new = c('eid', 'Age_AC', 'Smoking', 'Drinking'))
setnames(requirement, old = colnames(requirement), new = c('EID', 'Qualifications', 'Diabetes'))

# liver dysfunction = ALAT < 35 IU/l (male), < 25 IU/l (female)  | ASAT < 35 IU/l (male), < 25 IU/l (female)  |	GGT < 50 IU/l (male), < 35 IU/l (female)
# Modify row values for specific columns :

baseline_df$Smoking <- str_replace_all(baseline_df$Smoking, c("-3" = "Preferd no answer", "2" = "Current","0" = "No", "1" = "Yes" ))
baseline_df$Drinking <- str_replace_all(baseline_df$Drinking, c("-3" = "Preferd no answer", "2" = "Current", "0" = "No", "1" = "Yes" ))
baseline_df$Diabetes <- str_replace_all(baseline_df$Diabetes, c("-3" = "Preferd no answer", "-1" = "Do not know", "0" = "No", "1" = "Yes" ))
baseline_df$Qualifications <- str_replace_all(baseline_df$Qualifications, c("1" = "College/University", "2" = "A/AS-levels", "3" = "O-levels/GCSE", "4" = "CSE", "5" = "NVQ/HND/HNC", "6" = "Other-eg:nursing/teaching", "-7" = "None of above", "-3" = "Prefer not to answer"))

# Filter specific values for columns :

requirement$Diabetes <- ifelse(requirement$Diabetes == 0 | requirement$Diabetes == 1, requirement$Diabetes, NA)
covariates2$Drinking <- ifelse(covariates2$Drinking == 0 | covariates2$Drinking == 1, covariates2$Drinking, NA)
covariates2$Smoking <- ifelse(covariates2$Smoking == 0 | covariates2$Smoking == 1, covariates2$Smoking, NA)

# Remove rows with NA values :

requirement <- na.omit(requirement)
covariates2 <- na.omit(covariates2)

# Replace specific values with NA in data frame columns : 

for (col in c("EID", "Qualifications", "Diabetes")) { bp1[bp1[[col]] %in% c("Preferd no answer", "Do not know"), col] <- NA}
for (col in c("eid", "Age_AC", "Smoking", "Drinking")) { bp1[bp1[[col]] %in% c("Preferd no answer", "Do not know", "never", "Prefer not to answer", "Current"), col] <- NA }

# Merge data frames based on common IDs :

common7 <- intersect(covariates2$eid, requirement$EID)
df7_common <- covariates2[common7, ]
diabetes_common <- requirement[common7, ]
baseline_p <- merge(requirement, covariates2, by.x = 'EID', by.y = 'eid')

common8 <- intersect(baseline_p$EID, metabolites_baseline$Eid)
baseline_p_common <- baseline_p[common8, ]
metabolites_baseline_common <- metabolites_baseline[common8, ]
baseline <- merge(baseline_p, metabolites_baseline, by.x = 'EID', by.y = 'Eid')

bp1 <- intersect(baseline$EID, K760$eid_1)
baseline_common <- baseline[bp1, ]
K760_common <- K760[bp1, ]
bp1 <- merge(baseline, K760, by.x = 'EID', by.y = 'eid_1') %>%
  select("EID", "Qualifications", "Diabetes", "Age_AC", "Smoking", "Drinking", "BMI", "Sex", "Ethnicity", "Diagnosis") 

# Filter data based on specific values in a column:

filtered_bp1 <- baseline_df[baseline_df$Diagnosis == "K760", ]
filtered_bp2 <- baseline_df[baseline_df$Diagnosis != "K760", ]

selected_values <- c("K70.0", "K71.0", "K72.0", "K73.0", "K74.0", "K75.0", "K76.0", "K77.0", "C22.0","K75.8")
diagnosed_baseline <- baseline_df[baseline_df$Diagnosis_grp %in% selected_values, ]
healthy_baseline <- baseline_df[!(baseline_df$Diagnosis_grp %in% selected_values), ]

vtree(filtered_bp1, c("Diagnosis","Diabetes", "Drinking","Smoking"))
vtree(filtered_bp2, c("Diagnosis", "Diabetes", "Drinking", "Smoking"), maxNodes = 503000)

# Create CSV files: 

write.csv(baseline, "C:/Users/User/Desktop/PhD Documentation/My drafts/baseline_examination.csv", row.names = FALSE)
write.csv(bp2, "C:/Users/User/Desktop/Data/baseline_data.csv", row.names = FALSE)
write.csv(baseline_df, "C:/Users/User/Desktop/Data/Results/baseline.csv", row.names = FALSE)
write.csv(covariates, "C:/Users/User/Desktop/Data/Results/covariates.csv", row.names = FALSE)
write.csv(nutritrion, "C:/Users/User/Desktop/Data/nutritrion.csv", row.names = FALSE)
write.csv(init_exam, "C:/Users/User/Desktop/Data/init.csv", row.names = FALSE)