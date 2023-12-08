# Baseline files : 

baseline <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/baseline_data.csv')

# File Modification : 

setnames(baseline,old=colnames(baseline), new = c('Eid',"Education","Diabetes","Age","Smoking","Drinking","ICD_code"))

# Medication files : 

common9 <- intersect(baseline$Eid, medication$eid )
baseline_common = baseline[common9, ]
medication_common = medication[common9, ]
df9 = merge(baseline,medication, by.x=c('Eid'), by.y=c('eid'))

# Metformin 

metformin_cols <- grep("metformin", tolower(names(df9)), value = TRUE)
metformin_prop <- df9[, c('Eid',"Education","Diabetes","Age","Smoking","Drinking","ICD_code","Medications_1","Medications_2","Medications_3","Medications_4","Medications_5","Medications_6","Medications_7","Medications_8","Medications_9","Medications_10","Medications_11","Medications_12","Medications_13","Medications_14","Medications_15","Medications_16","Medications_17","Medications_18","Medications_19", metformin_cols)]
metformin_prop1 <- metformin_prop[, c("Eid", "Education", "Diabetes", "Age", "Smoking", "Drinking", "ICD_code","Medications_1","Medications_2","Medications_3","Medications_4","Medications_5","Medications_6","Medications_7","Medications_8","Medications_9","Medications_10","Medications_11","Medications_12","Medications_13","Medications_14","Medications_15","Medications_16","Medications_17","Medications_18","Medications_19", metformin_cols)]

# NASH with Metformin medicine

metformin_NASH <- metformin_prop[metformin_prop$ICD_code == "K758", ]

# MAFLD with Metformin medicine

metformin_MAFLD <- metformin_prop[metformin_prop$ICD_code == "K760", ]

# Non-metformin

non_metformin_cols <- grepl("glitazone|Pioglitazone|Liraglutide|Semaglutide|SGLT2 Inhibitors|Rapamycin|Clopidogrel|Beta blocker|Fenoterol|Vitamin D|ramipril", tolower(names(df9)))
non_metformin_cols <- names(df9)[non_metformin_cols]
non_metformin_prop <- df9[, c('Eid',"Education","Diabetes","Age","Smoking","Drinking","ICD_code","Medications_1","Medications_2","Medications_3","Medications_4","Medications_5","Medications_6","Medications_7","Medications_8","Medications_9","Medications_10","Medications_11","Medications_12","Medications_13","Medications_14","Medications_15","Medications_16","Medications_17","Medications_18","Medications_19",non_metformin_cols)]

# MAFLD without Metformin medicine 

non_metformin_NASH <- metformin_prop[non_metformin_prop$ICD_code == "K758", ]

# MAFLD without Metformin medicine

non_metformin_MAFLD <- metformin_prop[non_metformin_prop$ICD_code == "K760", ]

#File generation

write.csv(metformin_NASH, "C:/Users/User/Desktop/Data/Results/non_metformin_NASH.csv", row.names=FALSE)
write.csv(metformin_MAFLD, "C:/Users/User/Desktop/Data/Results/non_metformin_MAFLD.csv", row.names=FALSE)