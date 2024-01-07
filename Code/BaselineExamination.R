# Death data with baseline examination: 

# Read the required data files :

requirement <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/baselineRequire.csv') %>%
  select('eid', 'X6138.0.0', 'X2443.0.0')
covariates2 <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates2.csv') %>%
  select('eid', 'X21003.0.0', 'X20116.0.0', 'X20117.0.0')
metabolites_baseline <- metabolites %>%
  select('Eid', 'BMI', 'Sex', 'Ethnicity')

# Read baseline data file :

baseline_df <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/baseline_data.csv')

# Rename column names :

setnames(baseline_df, old = colnames(baseline_df), new = c('eid', 'Qualifications', 'Diabetes', 'Age_AC', 'Smoking', 'Drinking', 'Diagnosis'))
baseline_df<- as.data.frame(baseline_df[!duplicated(baseline_df$eid),])

mafld_baseline <- baseline_df[Diagnosis == 'K760', ]
mafld_baseline <- baseline_df[baseline_df$Diagnosis == 'K760', ]
nash_baseline <- baseline_df[baseline_df$Diagnosis == 'K758', ]
non_mafld_baseline <- baseline_df[baseline_df$Diagnosis != 'K760', ]
non_nash_baseline <- baseline_df[baseline_df$Diagnosis != 'K758', ]

# Read another data file :

baseline_df <- fread("C:/Users/User/Desktop/Data/ukb52200.csv") %>%
  select("eid", "X21003.0.0", "X20116.0.0", "X20117.0.0", "X41270.0.0", "X6138.0.0", "X2443.0.0")

hypertension_baseline_df <- fread("C:/Users/User/Desktop/Data/ukb52200.csv") %>%
  select("eid", "21003-0.0", "20116-0.0", "20117-0.0", "41270-0.0", "6138-0.0", "2443-0.0", "4080-0.0","4079-0.0")
Covariates2<- fread("ukb52200.csv", select=c("eid", "21000-0.0","21001-0.0"))

# Rename column names :

setnames(covariates2, old = colnames(covariates2), new = c('eid', 'Age_AC', 'Smoking', 'Drinking'))
setnames(requirement, old = colnames(requirement), new = c('EID', 'Qualifications', 'Diabetes'))

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

for (col in c("EID", "Qualifications", "Diabetes")) {
  bp1[bp1[[col]] %in% c("Preferd no answer", "Do not know"), col] <- NA
}

for (col in c("eid", "Age_AC", "Smoking", "Drinking")) {
  bp1[bp1[[col]] %in% c("Preferd no answer", "Do not know", "never", "Prefer not to answer", "Current"), col] <- NA
}

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