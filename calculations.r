# Dataframe 

df9 = df5%>%
  select('Eid_1','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19')
medication_group = df3%>%
  select('Eid','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19')
df7 = df5[df5$icd_group =='K71.0',] %>%
  select('Age','Sex','BMI','Ethnicity','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19')
medication_group = medication_group[!duplicated(medication_group), ]

# BMI as Patient Parameter

df7[df7$BMI > 18.5 & df7$BMI <= 25, "BMI_group"] <- "18.5-25"
df7[df7$BMI > 25 & df7$BMI <= 30, "BMI_group"] <- "25-30"
df7[df7$BMI > 30, "BMI_group"] <- "> 30"

# Age as Patient Parameter

df7[df7$Age <= 12, "age_group"] <- "0-12"
df7[df7$Age > 13 & df7$Age <= 19, "age_group"] <- "13-19"
df7[df7$Age > 20 & df7$Age <= 35, "age_group"] <- "20-35"
df7[df7$Age > 35 & df7$Age <= 50, "age_group"] <- "35-50"
df7[df7$Age > 50 & df7$Age <= 65, "age_group"] <- "50-65"
df7[df7$Age > 65, "age_group"] <- "> 65"

# Healthy Patient criteria

Hepa = df5%>%
  select('Eid_1','diag_icd10','Total_Free_C','Liver_fat','HDL_C','Total_Triglycerides','Clinical_LDL_C','icd_group')

# baseline metabolites enzymes parameters for new dataframe Healthy Patients. 

max(Hepa$Concentration_of_HDL, na.rm = TRUE)
Hepa1 = filter(Hepa, Total_Triglycerides < 1.5, )
Hepa1 = Hepa1[!duplicated(Hepa1), ]
Hepa2 = filter(Hepa, Total_Triglycerides < 1.5 & Total_Free_C < 2, )
Hepa2 = Hepa2[!duplicated(Hepa2), ]
Hepa3 = filter(Hepa, Total_Triglycerides < 1.5 & Total_Free_C < 2 & Liver_fat < 5, )
Hepa3 = Hepa3[!duplicated(Hepa3), ]
Hepa_T = subset(Hepa,Liver_fat < 5 & HDL_C < 1.55 & Clinical_LDL_C < 2.6 & Total_Triglycerides < 1.5 & Total_Free_C < 2, )
#& ALAT < 35 & ASAT < 35, GGT < 50



# Extracting Disease specific codes

Hepa_T <- subset(Hepa, !(icd_group %in% c("K70.0", "K71.0", "K72.0", "K73.0", "K74.0", "K75.0", "K76.0", "K77.0","C22.0")))
Hepa_T = Hepa_T[!duplicated(Hepa_T), ]


# Mean values from Healthy Patients df: 

Hepa_T_mean <- colMeans(Hepa_T[, c(3, 4, 5, 6, 7, 8, 9)])

# 1. Standard deviation Min and Max Values
# calculate summary statistics for each column

metabolites2 <- metabolites %>%
  select('Total_C','Total_C-HDL_C','Remnant_C','VLDL_C','Clinical_LDL_C','LDL_C','HDL_C','Total_Triglycerides','Triglycerides_in_VLDL', 'Triglycerides_in_LDL','Triglycerides_in_HDL','Total_Phospholipids_in_Lipoprotein','Phospholipids_in_VLDL', 'Phospholipids_in_LDL','Phospholipids_in_HDL','Total_Esterified_C', 'Cholesteryl_Esters_in_VLDL','Cholesteryl_Esters_in_LDL','Cholesteryl_Esters_in_HDL','Total_Free_C','FC_in_VLDL','FC_in_LDL','FC_in HDL', 'Total Lipids in Lipoprotein ', 'Total Lipids in VLDL','Total Lipids in LDL','Total Lipids in HDL','Total  Lipoprotein ','VLDL','LDL','Concentration_of_HDL','Average Diameter for VLDL ','Average Diameter for LDL ,Average Diameter for HDL_','Phosphoglycerides','Triglycerides_to_Phosphoglycerides ratio','Total_Cholines','Phosphatidylcholines','Sphingomyelins','Apolipoprotein_B','Apolipoprotein_A1','Apolipo protein B-Apolipoprotein A1 ratio',
         'Total_Fatty_Acids','Degree of Unsaturation','Omega-3 Fatty_Acids','Omega-6_Fatty_Acids','Polyunsaturated_Fatty_Acids','Monounsaturated_Fatty_Acids','Saturated_Fatty_Acids','Linoleic_Acid','DHA','Omega-3_Fatty_Acids_to_Total_Fatty_Acids_percentage','Omega-6 Fatty Acids to Total Fatty Acids percentage','Polyunsaturated Fatty Acids to Total Fatty Acids percentage','Monounsaturated Fatty Acids to Total Fatty Acids percentage','Saturated Fatty Acids to Total Fatty Acids percentage','Linoleic Acid to Total Fatty Acids percentage','Docosahexaenoic Acid to Total Fatty Acids percentage','Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio','Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio','Alanine','Glutamine','Glycine','Histidine','Total  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 Chylomicrons and EL_VLDL ','Total Lipids in Chylomicrons and EL_VLDL',
         'Phospholipids in Chylomicrons and EL_VLDL','C in Chylomicrons and EL_VLDL','CE_ in Chylomicrons and EL_VLDL','FC_in Chylomicrons and EL_VLDL','Triglycerides in Chylomicrons and EL_VLDL',' VL_VLDL Particle','Total Lipids in VL_VLDL','Phospholipids in VL_VLD','C_ in VL_VLDL','CE_ in VL_VLDL','	FC_in VL_VLDL','Triglycerides in VL_VLDL',' L_VLDL ','Total Lipids in L_VLDL','Phospholipids in L_VLDL','C_ in L_VLDL','CE_ in L_VLDL','FC_in L_VLDL','Triglycerides in L_VLDL','M_VLDL','Total Lipids in M_VLDL','Phospholipids in M_VLDL','C_ in M_VLDL','CE_ in M_VLDL','FC_in M_VLDL','Triglycerides in M_VLDL','S_VLDL','Total Lipids in S_VLDL','Phospholipids in S_VLDL','C_ in S_VLDL','CE_ in S_VLDL','FC_in S_VLDL','Triglycerides in S_VLDL',
         'VS_VLDL','Total Lipids in VS_VLDL','Phospholipids in_VS_VLDL','C_ in_VS_VLDL','CE_ in VS_VLDL','FC_in VS_VLDL','Triglycerides in VS_VLDL','IDL','Total Lipids in IDL','Phospholipids in IDL','C_ in IDL','	CE_ in IDL','FC_in_IDLs','Triglycerides in IDL','L_LDL','	Total Lipids in L_LDL','Phospholipids in L_LDL','C_ in L_LDL','CE_ in L_LDL','FC_in L_LDL','Triglycerides in L_LDL',' M_LDL ','Total Lipids in M_LDL','Phospholipids in M_LDL','C_ in M_LDL','CE_ in M_LDL','FC_in M_LDL','Triglycerides in M_LDL',' S_LDL ','Total Lipids in S_LDL','Phospholipids in S_LDL','C_ in S_LDL','CE_ in S_LDL','FC_in S_LDL','Triglycerides in S_LDL',' VL_HDL','Total Lipids in VL_HDL','Phospholipids in VL_HDL','C_ in VL_HDL',
         'CE_ in VL_HDL','FC_in VL_HDL','Triglycerides in_VL_HDL',' L_HDL_','Total Lipids in L_HDL','Phospholipids in L_HDL','C_ in L_HDL','CE_ in L_HDL','FC_ in L_HDL','Triglycerides in L_HDL',' M_HDL','Total Lipids in M_HDL','Phospholipids in M_HDL','C_ in M_HDL')
metabolites2 <- data.frame(sapply(metabolites2, as.numeric))

metabolites2[metabolites2 == -999] <- NA

# Create summary table for standard deviation, Internal quartile range with(0.25,0.75)

HepaT = Hepa_T%>%
  select('Eid_1','diag_icd10','Total_Free_C','Liver_fat','HDL_C','Total_Triglycerides','Clinical_LDL_C')

# 'ALAT','ASAT','GGT' genes for observing summary.  

summary_table <- data.frame(
  Variables = names(metabolites),
  Mean = round(colMeans(metabolites, na.rm = TRUE), 2),
  SD = round(sapply(metabolites, sd, na.rm = TRUE), 2),
  Median = round(sapply(metabolites, median, na.rm = TRUE), 2),
  Q1 = round(sapply(metabolites, quantile, probs = 0.25, na.rm = TRUE), 2),
  Q3 = round(sapply(metabolites, quantile, probs = 0.75, na.rm = TRUE), 2),
  IQR = round(sapply(metabolites, IQR, na.rm = TRUE), 2),
  Min = round(sapply(metabolites, min, na.rm = TRUE), 2),
  Max = round(sapply(metabolites, max, na.rm = TRUE), 2),
  stringsAsFactors = FALSE)

# Create a function to calculate the summary statistics

Hepa_TT <- Hepa_TT[, -c(which(names(Hepa_TT) == "Eid_1"), which(names(Hepa_TT) == "diag_icd10"), which(names(Hepa_TT) == "icd_group"))]
Hepa_TT <- as.data.frame(lapply(Hepa_T, as.numeric))

# Healthy patient table to analyze for mean, median, IQR dataframe. 

Hepa_TT_summary <- data.frame(
  Variables = names(Hepa_TT),
  sd = round(apply(Hepa_TT, 2, sd, na.rm = TRUE), 2),
  IQR = round(apply(Hepa_TT, 2, IQR, na.rm = TRUE), 2),
  min = round(apply(Hepa_TT, 2, min, na.rm = TRUE), 2),
  Q1 = round(apply(Hepa_TT, 2, quantile, probs = 0.25, na.rm = TRUE), 2),
  median = round(apply(Hepa_TT, 2, median, na.rm = TRUE), 2),
  Q3 = round(apply(Hepa_TT, 2, quantile, probs = 0.75, na.rm = TRUE), 2),
  max = round(apply(Hepa_TT, 2, max, na.rm = TRUE), 2),
  mean = round(apply(Hepa_TT, 2, mean, na.rm = TRUE), 2))

# Non Medications Group for healthy patients dataset. 

No_Medication = filter(medication_group, rowSums(is.na(medication_group)) == ncol(medication_group)-1)
Healthy_Patients <- No_Medication[!duplicated(No_Medication), ]

# Write created dataframes for saving generated content 

write.csv(No_medication_1, "C:/Users/User/Desktop/PhD Documentation/My drafts/No_Medication.csv", row.names=FALSE)
write.csv(Hepa_T, "C:/Users/User/Desktop/PhD Documentation/My drafts/Healthy_Patients.csv", row.names=FALSE)
write.csv(summary_table, "C:/Users/User/Desktop/PhD Documentation/My drafts/summary_table.csv", row.names=FALSE)
write.csv(death_info, "C:/Users/User/Desktop/PhD Documentation/My drafts/death_info.csv", row.names=FALSE)

#Graphical Observaiention

vtree(df5, c("diag_icd10","Sex","BMI_group","age_group"))
vtree(df7, c("diag_icd10","Sex","BMI_group","age_group"))
pie(occurance2, by.x=c(Freq), by.y=c(Var1))