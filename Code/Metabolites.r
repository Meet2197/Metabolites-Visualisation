# Libraries to upload file for data modifications, generating data table, plotting the data. 

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(data.table)
library(readr)
library(janitor)
library(DT)
library(naniar)
library(ggplot2)

# Read Text and CSV Files_the R Script 

medication <- read.csv(file = 'C:/Users/User/Desktop/Data/medication.csv')[ ,2:21]
covariates <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates_Markus.csv')[ ,2:7]
hesin_diag <- read.table(file = 'C:/Users/User/Desktop/Data/hesin_diag.txt' ,header = TRUE , sep = "\t" )%>% 
  select('requires.eid','ins_index','diag_icd10')
hesin <- read.table(file = 'C:/Users/User/Desktop/Data/hesin.txt', header = TRUE , sep = "\t") %>% 
  select('eid','ins_index','epistart','epiend','disdate','admidate')
metabolites <- read.csv(file = 'C:/Users/User/Desktop/Data/metabolites.csv')
mri <- read.csv(file = 'C:/Users/User/Desktop/Data/MRI.csv')%>%
  select('eid','X40061.2.0')
death <- read.table(file = 'C:/Users/User/Desktop/Data/Endpoints/death/death.txt', header = TRUE , sep = "\t")%>% 
  select('eid','ins_index','date_of_death')
death_cause <- read.table(file = 'C:/Users/User/Desktop/Data/Endpoints/death/death_cause.txt', header = TRUE , sep = "\t")%>% 
  select('eid','ins_index','cause_icd10')
covariates2 <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates_Meet2.csv') %>%
  select('eid','X30620.0.0','X30650.0.0','X30730.0.0')
genes <- read.csv(file = 'C:/Users/User/Desktop/Data/genes.csv') %>%
  select("eid","rs2642438_A","rs1260326_T","rs780094_T","rs72613567_TA","rs10069690_T","rs9272105_A",
         "rs2856718_T","rs7453920_A","rs3077_G","rs9277535_G","rs58542926_T","rs455804_A","rs738409_G")   
NASH <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/Nash.csv')

# After setnames of metabolites,

metabolite <- metabolites %>%
  select('Eid','Ethnicity','BMI','Age','Sex')

# 'Total_Free_C','Clinical_LDL_C','HDL_C','Total_Triglycerides','Linoleic_Acid','Omega-3 _FA_','Omega-6__FA_','Triglycerides_VLDL', 'Triglycerides_LDL','Triglycerides_HDL','CE_CM and EL_VLDL','CE_VL_VLDL','CE_L_VLDL','CE_M_VLDL','CE_S_VLDL','CE_VS_VLDL','	CE_IDL','CE_L_LDL','CE_M_LDL','CE_S_LDL','CE_VL_HDL','CE_L_HDL','CE_L_HDL',

# Mutation of data for Medication
 
setnames(hesin_diag,"requires.eid",'eid_1')
setnames(covariates, old=colnames(covariates), new = c('eid','Birth_Year','Birth_Month','BMI','Age_at_Obs.','Gender'))
setnames(medication, old=colnames(medication), new = c('eid','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19'))
setnames(metabolites, old=colnames(metabolites), new = c('No.','Eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Concentration_of_HDL','Average Diameter for VLDL ','Average Diameter for LDL ,Average Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
                                                         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
                                                         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
                                                         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
                                                         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','Ethnicity','BMI','Age','Sex'))
setnames(mri,old=colnames(mri), new = c('eid',"Liver_fat"))
setnames(covariates2,old=colnames(covariates2),new = c('eid','ALAT','ASAT','GGT'))
setnames(death,"eid","eid_1")
setnames(NASH,"X41202.0.0","icd_code")
setnames(K760,"K760","Diagnosis")


# sapply(data2$epistart, typeof), class(data2$epistart)

hesin$epistart <- dmy(hesin$epistart)
hesin$epiend <- dmy(hesin$epiend)
hesin$admidate <- dmy(hesin$admidate)
hesin$disdate <- dmy(hesin$disdate)

# Days difference

hesin$epi_diff = difftime(hesin$epiend,hesin$epistart,units = "days")
hesin$admi_diff = difftime(hesin$disdate,hesin$admidate,units = "days")

# creating round of values of metabolites and diabetes. 

metabolites <- round(na.omit(metabolites) , 3)
covariates <- na.omit (covariates)
mri <- na.omit(mri)
diabetes <- na.omit(diabetes)

# Sellectign columns form admissionP dataframe. 

admissionP <- na.omit(hesin) %>%
  select("eid","admidate")

# death info consists several columns mentioned below for MAFLD and NASH associated death.  

death_info <- na.omit(df6) %>%
  select("eid_1","date_of_death", "icd_group","icd_names")

colnames(death_info)

medication[is.na(medication)] <- simvastatin
medication[medication == simvastatin] <- 0
medication[medication == simvastatin] <- NA

# 1.1 Find Common Ids_the Data frame:

common <- intersect(hesin$eid, hesin_diag$eid_1)
hesin_common = hesin[common, ]
hesin_diag_common = hesin_diag[common, ]
df1 = merge(hesin_diag,hesin, by.x=c('eid_1','ins_index'), by.y=c('eid','ins_index'))%>% 
  select('eid_1','diag_icd10', 'admi_diff', 'epi_diff', 'icd_group')

common2 <- intersect(df1$eid_1, metabolites$Eid)
df2_common = df1[common2, ]
metabolites_common = metabolites[common2, ]
df2 = merge(df1, metabolites, by.x=c('eid_1'), by.y=c('Eid')) %>%
  select('eid_1','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL','diag_icd10','icd_group')

## 1.2 Medications dataframe generation. 

common3 <- intersect(df2$eid_1, medication$eid)
df3_common = df2[common3, ]
medication_common = medication[common3, ]
df3 = merge(medication, df2, by.x=c('eid'), by.y=c('eid_1')) %>%
  select('eid','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19','diag_icd10','icd_group')

df6 = merge(medication, df2, by.x=c('eid'), by.y=c('eid_1')) %>%
  select('eid','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL','diag_icd10','icd_group')
         
# NASH icd code comaprison with Sub metabolites. 

Nash = merge(NASH, metabolites, by.x=c('eid'), by.y=c('Eid')) %>%
  select('eid','icd_code','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL')

common5 <- intersect(df2$eid_1, K758$eid_1)
df5_common = df2[common5, ]
K758_common = K758[common5, ]
df5 = merge(df2, K758, by.x=c('eid_1'), by.y=c('eid_1')) %>%
  select('eid_1','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL','diag_icd10') 

MAFLD_df <- subset(df2, startsWith(as.character(diag_icd10), 'K760'))
Fibrosis_cirrhosis <- df2[df2$icd_group == 'K74.0' & !is.na(df2$icd_group), ]
neoplasm <- df2[df2$icd_group == 'C22.0' & !is.na(df2$icd_group), ]
no_liver_disease <- df2[!(df2$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0')), ]

# Nash (K75.8) dataframe with metabolites columns to observe impact over NASH development.  

NASH_dataframe <- intersect(NASH$eid, metabolites$Eid)
nash_common = NASH[NASH_dataframe, ]
metabolites_common = metabolites[NASH_dataframe, ]
NASH_dataframe = merge(NASH,metabolites, by.x=c('eid'), by.y=c('Eid'))%>% 
  select('eid','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','P_HDL')
NASH_df = NASH_dataframe[!duplicated(NASH_dataframe), ]

# No medication Dataframe generated_calculations.R file with no usage of medications for NASH and MAFLD patients. 

common4 <- intersect(df2$eid_1, No_Medication$eid)
df4_common = df2[common4, ]
No_Medication_common = No_Medication[common3, ]
healthy_df = merge(df2, No_Medication, by.x=c('eid_1'), by.y=c('eid')) %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL')

healthy_df = healthy_df[!duplicated(healthy_df), ]
count(metformin_df)

write.csv(metformin_table, "C:/Users/User/Desktop/Data/Results/metformin_table.csv", row.names=FALSE)
write.csv(negative_control_data_table, "C:/Users/User/Desktop/Data/Results/non_metformin_dt.csv", row.names=FALSE)


healthy_metformin_df <- metformin_df[!(metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0')), ]
diagnosed_metformin_df <- metformin_df[metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0'), ]
MAFLD_df <- diagnosed_metformin_df[diagnosed_metformin_df$icd_group == 'K76.0', ]

#Frequency
frequency_MAFLD_nonmetformin_df <- (nrow(MAFLD_df) / count) * 100
NASH_df <- metformin_df[metformin_df$diag_icd10 == 'K758', ]
frequency_NASH_df <- (nrow(MAFLD_df) / 91338) * 100 

# non-metformin : 

non_metformin_rows <- df3[rowSums(df3[, grepl("Medications_", names(df3))] == "metformin", na.rm = TRUE) == 0, ]
non_metformin_df <- non_metformin_rows[, c('eid', 'icd_group', 'diag_icd10', grep("Medications_", names(df3), value = TRUE))]
non_metformin_df = non_metformin_df[!duplicated(non_metformin_df), ]
healthy_nonmetformin_df <- non_metformin_df[!(non_metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0')), ]
diagnosed_nonmetformin_df  <- non_metformin_df[non_metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0'), ]
MAFLD_nonmetformin_df <- diagnosed_nonmetformin_df[diagnosed_nonmetformin_df$icd_group == 'K76.0', ]
frequency_MAFLD_nonmetformin_df <- (nrow(MAFLD_nonmetformin_df) / 91338) * 100

## 1.3 Covariates file include baseline with intersect with df3. 

common4 <- intersect(df3$Eid, covariates$Eid_1)
df3_common = df3[common4, ]
covariates_common = covariates[common4, ]
df4 = merge(covariates, df3, by.x=c('Eid_1'), by.y=c('Eid'))

## 1.4 MRI data intersect df3.

common5 <- intersect(df4$Eid_1 , mri$eid)
df4_common = df4[common5, ]
mri_common = mri[common5,]
df5 = merge(df3, mri, by.x=c('Eid_1'), by.y=c('eid'))

## 1.5 Death cause and death dataframe intersect to observe cause of NAFLD and MAFLD. 

common6 <- intersect(death$eid_1 , death_cause$eid)
death_common = death[common6, ]
death_cause_common = death_cause[common6,]
death_df = merge(death, death_cause, by.x=c('eid_1','ins_index'), by.y=c('eid','ins_index'))
NASH_death <- death_df[death_df$cause_icd10 == "K758", ]
MASLD_death <- death_df[death_df$cause_icd10 == "K760", ]

# Filter df6 based on icd_group

MAFLD_Death_filtered <- df6[df6$icd_group %in% c("K76.0"), ]
NASH_Death_filtered <- df6[df6$cause_icd10 == "K758", ]

# Display the results
print(MAFLD_Death_filtered)

## 1.7 Genes file intersect with df6. 

common7 <- intersect(df6$eid_1, genes$eid)
death_genetic_common = df6[common7, ]
death_cause_common = genes[common7,]
df7 = merge(df6, genes, by.x=c('eid_1'), by.y=c('eid'))

common8 <- intersect(df7$eid_1, metabolite$Eid)
death_genetic_common = df7[common8, ]
metabolite_common = metabolite[common8,]
df8 = merge(df7, metabolite, by.x=c('eid_1'), by.y=c('Eid'))

genetic_compare <- c("K700", "K710", "K720", "K730", "K740", "K750", "K760", "K770","C220")
genetic_compare_df <- df6[df6$date_of_death %in% genetic_compare, ]

# Medications frequent occurrence_Medications  

occurance <- data.frame(table(medication$Medications_2))
occurance2 = occurance[occurance$Freq > 1000,]

# Use mutate_all and recode to transform all row values

# for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140875408] <- "allopurinol" } 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140883504] <- "cetirizine" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140916682] <- "evening primrose oil"} 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140910814] <- "sodium thyroxine" } 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140863152] <- "diazepam" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140865580] <- "asacol-400mg-e/c-tablet" } 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871310] <- "ibuprofen" } 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140868226] <- "Aspirin" } 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879616] <- "amitriptyline" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140922174] <- "alendronate_sodium"} 
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140870390] <- "ferrous sulphate" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871310] <- "ibuprofen" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140888366] <- "thiamine_preparation" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140860696] <- "lisinopril" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140864752] <- "lansoprazole" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871024] <- "vitamin_b_compound_tablet" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879760] <- "bisoprolol" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879802] <- "amlodipine" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140923346] <- "co-codamol" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141182628] <- "tiotropium" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141176832] <- "seretide_50_evohaler" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141145660] <- "valsartan" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140863144] <- "zopiclone" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140926606] <- "salbutamol_100micrograms_spacehaler" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140888510] <- "verapamil" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140883066] <- "insulin product" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140860806] <- "ramipril" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140909708] <- "furosemide" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 99999] <- "uncoded" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140864992] <- "tramadol" }
for (i_seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140884600] <- "metformin" }
for (i in seq_along(medication)) {
  medication[[i]][medication[[i]] == 1140884600] <- "metformin"
}
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140865716] <- "senna" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 2038460150] <- "paracetamol" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140874930] <- "prednisolone" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140888266] <- "warfarin" } 
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140910766] <- "nicorandil" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140874420] <- "quinine" } 
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140879406] <- "ranitidine" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140881856] <- "salbutamol" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141191044] <- "levothyroxine sodium" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 6] <- "Non-cancer Illness"}
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 3004] <- "Mobility problem severity" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140851028] <- "chalk" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140884412] <- "sumatriptan" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140865634] <- "omeprazole" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 2038460150] <- "paracetamol" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140866738] <- "atenolol" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141146234] <- "NA" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140861958] <- "NA" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141194794] <- "bendroflumethiazide" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141188442] <- "glucosamine product" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140909674] <- "cod_liver_oil_capsule" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140868226] <- "Aspirin" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140852976] <- "multivitamins" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140861998] <- "ventolin 100micrograms inhaler" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1193] <- "omega-3/fish oil supplement" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140909726] <- "vitamin c product" }
for (i_seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140888538] <- "zinc product" }

# Ethnicity values change 

for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1001] <- "British" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 6] <- "Other ethnic group" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3] <- "Asian" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1002] <- "Irish" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3001] <- "Indian" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% -3] <- "Prefer not to answer" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 4001] <- "Caribbean" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1003] <- "Any other white background" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1002] <- "Irish" }
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 2001] <- "White and Black Caribbean" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3002] <- "Pakistani" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 4002] <- "African" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3004] <- "Any other Asian background" } 
for (i_seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% -1] <- "Do not know" }

# Genes changes
for (i_seq_along(genes)) { genes[[i]][genes[[i]] %in% 0] <- "Physical activity" } %>% { genes[[i]][genes[[i]] %in% 1] <- "heterozygous" } %>% { genes[[i]][genes[[i]] %in% 2] <- "homozygous" }

#Hesin diagnosis coding

hesin_diag[hesin_diag$diag_icd10 %in% c("K700", "K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "icd_group"] <- "K70.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K710", "K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718", "K719"), "icd_group"] <- "K71.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K720", "K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728", "K729"), "icd_group"] <- "K72.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K730", "K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "icd_group"] <- "K73.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K740", "K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "icd_group"] <- "K74.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K750", "K751", "K752", "K753", "K754", "K755", "K756", "K757", "K759"), "icd_group"] <- "K75.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K760", "K761", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769"), "icd_group"] <- "K76.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K770", "K771", "K772", "K773", "K774", "K705", "K776", "K777", "K778", "K779"), "icd_group"] <- "K77.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K220", "K221", "K222", "K223", "K224", "K705", "K226", "K227", "K228", "K229"), "icd_group"] <- "C22.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K75.8"), "icd_group"] <- "K75.8"

# ICD CODE grouping

death_cause[death_cause$cause_icd10 %in% c("K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "icd_group"] <- "K70.0"
death_cause[death_cause$cause_icd10 %in% c("K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718","K719"), "icd_group"] <- "K71.0"
death_cause[death_cause$cause_icd10 %in% c("K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728","K729"), "icd_group"] <- "K72.0"
death_cause[death_cause$cause_icd10 %in% c("K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "icd_group"] <- "K73.0"
death_cause[death_cause$cause_icd10 %in% c("K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "icd_group"] <- "K74.0"
death_cause[death_cause$cause_icd10 %in% c("K751", "K752", "K753", "K754", "K755", "K756", "K757"), "icd_group"] <- "K75.0"
death_cause[death_cause$cause_icd10 %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767","K768", "K769"), "icd_group"] <- "K76.0"
death_cause[death_cause$cause_icd10 %in% c("K771", "K772", "K773", "K774", "K705", "K776", "K777","K778", "K779"), "icd_group"] <- "K77.0"
death_cause[death_cause$cause_icd10 %in% c("C221", "C222", "C223", "C224", "C225", "C226", "C227","C228", "C229"), "icd_group"] <- "C22.0"
death_cause[death_cause$cause_icd10 == "K75.8", "icd_group"] <- "K75.8"

death_cause[death_cause$icd_group %in% c("K70.0"), "icd_names"] <- "Alcoholic liver disease"
death_cause[death_cause$icd_group %in% c("K71.0"), "icd_names"] <- "Toxic liver disease "
death_cause[death_cause$icd_group %in% c("K72.0","K73.0"), "icd_names"] <- "Hepatic failure"
death_cause[death_cause$icd_group %in% c("K74.0"), "icd_names"] <- "Fibrosis and cirrhosis"
death_cause[death_cause$icd_group %in% c("K75.0"), "icd_names"] <- "Inflammatory liver diseases"
death_cause[death_cause$icd_group %in% c("K75.8"), "icd_names"] <- "NASH"
death_cause[death_cause$icd_group %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767"), "icd_names"] <- "Chronic Hepatitis MAFLD"
death_cause[death_cause$icd_group %in% c("C22.0"), "icd_names"] <- "liver Malignant neoplasm"


baseline_df[baseline_df$Diagnosis %in% c("K700", "K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "Diagnosis_grp"] <- "K70.0"
baseline_df[baseline_df$Diagnosis %in% c("K710", "K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718", "K719"), "Diagnosis_grp"] <- "K71.0"
baseline_df[baseline_df$Diagnosis %in% c("K720", "K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728", "K729"), "Diagnosis_grp"] <- "K72.0"
baseline_df[baseline_df$Diagnosis %in% c("K730", "K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "Diagnosis_grp"] <- "K73.0"
baseline_df[baseline_df$Diagnosis %in% c("K740", "K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "Diagnosis_grp"] <- "K74.0"
baseline_df[baseline_df$Diagnosis %in% c("K750", "K751", "K752", "K753", "K754", "K755", "K756", "K757", "K759"), "Diagnosis_grp"] <- "K75.0"
baseline_df[baseline_df$Diagnosis %in% c("K760", "K761", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769"), "Diagnosis_grp"] <- "K76.0"
baseline_df[baseline_df$Diagnosis %in% c("K770", "K771", "K772", "K773", "K774", "K705", "K776", "K777", "K778", "K779"), "Diagnosis_grp"] <- "K77.0"
baseline_df[baseline_df$Diagnosis %in% c("C220", "C221", "C222", "C223", "C224", "C705", "C226", "C227", "C228", "C229"), "Diagnosis_grp"] <- "C22.0"
baseline_df[baseline_df$Diagnosis == "K75.8", "Diagnosis_grp"] <- "K75.8"
