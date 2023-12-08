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

# Read Text and CSV Files in the R Script 

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
NASH <- read.csv(file = 'C:/Users/User/Desktop/Data/NASH.csv')

# After setnames of metabolites,

metabolite <- metabolites %>%
  select('Eid','Ethnicity','BMI','Age','Sex')

# 'Total_Free_C','Clinical_LDL_C','HDL_C','Total_Triglycerides','Linoleic_Acid','Omega-3 Fatty_Acids','Omega-6_Fatty_Acids','Triglycerides_in_VLDL', 'Triglycerides_in_LDL','Triglycerides_in_HDL','CE_ in Chylomicrons and EL_VLDL','CE_ in VL_VLDL','CE_ in L_VLDL','CE_ in M_VLDL','CE_ in S_VLDL','CE_ in VS_VLDL','	CE_ in IDL','CE_ in L_LDL','CE_ in M_LDL','CE_ in S_LDL','CE_ in VL_HDL','CE_ in L_HDL','CE_ in L_HDL',

# Mutation of data for Medication
 
setnames(hesin_diag,"requires.eid",'eid_1')
setnames(covariates, old=colnames(covariates), new = c('Eid_1','Birth_Year','Birth_Month','BMI','Age_at_Obs.','Gender'))
setnames(medication, old=colnames(medication), new = c('eid','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19'))
setnames(metabolites, old=colnames(metabolites), new = c('No.','Eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_in_VLDL', 'TG_in_LDL','TG_in_HDL','Tl_Phospholipids_in_Lipoprotein','Phospholipids_in_VLDL', 'Phospholipids_in_LDL','Phospholipids_in_HDL','Tl_Esterified_C', 'CE_in_VLDL','CE_in_LDL','CE_in_HDL','Tl_Free_C','FC_in_VLDL','FC_in_LDL','FC_in HDL', 'Tl Lipids in Lipoprotein ', 'Tl Lipids in VLDL','Tl Lipids in LDL','Tl Lipids in HDL','Tl  Lipoprotein ','VLDL','LDL','Concentration_of_HDL','Average Diameter for VLDL ','Average Diameter for LDL ,Average Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
                                                         'Tl_Fatty_Acids','Degree of Unsaturation','Omega-3 Fatty_Acids','Omega-6_Fatty_Acids','Polyunsaturated_Fatty_Acids','Monounsaturated_Fatty_Acids','Saturated_Fatty_Acids','Linoleic_Acid','DHA','Omega-3_Fatty_Acids_to_Tl_Fatty_Acids_percentage','Omega-6 Fatty Acids to Tl Fatty Acids percentage','Polyunsaturated Fatty Acids to Tl Fatty Acids percentage','Monounsaturated Fatty Acids to Tl Fatty Acids percentage','Saturated Fatty Acids to Tl Fatty Acids percentage','Linoleic Acid to Tl Fatty Acids percentage','Docosahexaenoic Acid to Tl Fatty Acids percentage','Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio','Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 Chylomicrons and EL_VLDL ','Tl Lipids in Chylomicrons and EL_VLDL',
                                                         'Phospholipids in Chylomicrons and EL_VLDL','C in Chylomicrons and EL_VLDL','CE_ in Chylomicrons and EL_VLDL','FC_in Chylomicrons and EL_VLDL','TG in Chylomicrons and EL_VLDL',' VL_VLDL Particle','Tl Lipids in VL_VLDL','Phospholipids in VL_VLD','C_in_VL_VLDL','CE_ in VL_VLDL','	FC_in VL_VLDL','TG in VL_VLDL','L_VLDL','Tl Lipids in L_VLDL','Phospholipids in L_VLDL','C_in_L_VLDL','CE_ in L_VLDL','FC_in L_VLDL','TG in L_VLDL','M_VLDL','Tl Lipids in M_VLDL','Phospholipids in M_VLDL','C_in_M_VLDL','CE_ in M_VLDL','FC_in M_VLDL','TG in M_VLDL','S_VLDL','Tl Lipids in S_VLDL','Phospholipids in S_VLDL','C_in_S_VLDL','CE_ in S_VLDL','FC_in S_VLDL','TG in S_VLDL',
                                                         'VS_VLDL','Tl Lipids in VS_VLDL','Phospholipids in_VS_VLDL','C_in_VS_VLDL','CE_ in VS_VLDL','FC_in VS_VLDL','TG in VS_VLDL','IDL','Tl Lipids in IDL','Phospholipids in IDL','C_in_IDL','CE_in_IDL','FC_in_IDLs','TG_in_IDL','L_LDL','Tl Lipids in L_LDL','Phospholipids in L_LDL','C_in_L_LDL','CE_in_L_LDL','FC_in L_LDL','TG in L_LDL','M_LDL','Tl Lipids in M_LDL','Phospholipids in M_LDL','C_ in M_LDL','CE_ in M_LDL','FC_in M_LDL','TG in M_LDL','S_LDL','Tl Lipids in S_LDL','Phospholipids in S_LDL','C_in_S_LDL','CE_ in S_LDL','FC_in S_LDL','TG in S_LDL',' VL_HDL','Tl Lipids in VL_HDL','Phospholipids in VL_HDL','C_in_VL_HDL',
                                                         'CE_in_VL_HDL','FC_in VL_HDL','TG_in_VL_HDL','L_HDL','Tl Lipids in L_HDL','Phospholipids in L_HDL','C_in_L_HDL','CE_in_L_HDL','FC_ in L_HDL','TG in L_HDL','M_HDL','Tl Lipids in M_HDL','Phospholipids in M_HDL','C_in_M_HDL','Ethnicity','BMI','Age','Sex'))
setnames(mri,old=colnames(mri), new = c('eid',"Liver_fat"))
setnames(covariates2,old=colnames(covariates2),new = c('eid','ALAT','ASAT','GGT'))
setnames(death,"eid","eid_1")
setnames(NASH,"X41202.0.0","icd_code")
setnames(K760,"K760","Diagnosis")


#sapply(data2$epistart, typeof), class(data2$epistart)

hesin$epistart <- dmy(hesin$epistart)
hesin$epiend <- dmy(hesin$epiend)
hesin$admidate <- dmy(hesin$admidate)
hesin$disdate <- dmy(hesin$disdate)

# Days difference

hesin$epi_diff = difftime(hesin$epiend,hesin$epistart,units = "days")
hesin$admi_diff = difftime(hesin$disdate,hesin$admidate,units = "days")

metabolites <- round(na.omit(metabolites) , 3)
covariates <- na.omit (covariates)
mri <- na.omit(mri)
diabetes <- na.omit(diabetes)

admissionP <- na.omit(hesin) %>%
  select("eid","admidate")

death_info <- na.omit(df6) %>%
  select("eid_1","date_of_death", "icd_group","icd_names")

colnames(death_info)

medication[is.na(medication)] <- simvastatin
medication[medication == simvastatin] <- 0
medication[medication == simvastatin] <- NA

# 1.1 Find Common Ids in the Data frame:

common <- intersect(hesin$eid, hesin_diag$eid_1)
hesin_common = hesin[common, ]
hesin_diag_common = hesin_diag[common, ]
df1 = merge(hesin_diag,hesin, by.x=c('eid_1','ins_index'), by.y=c('eid','ins_index'))%>% 
  select('eid_1','diag_icd10', 'admi_diff', 'epi_diff', 'icd_group')

common2 <- intersect(df1$eid_1, metabolites$Eid)
df2_common = df1[common2, ]
metabolites_common = metabolites[common2, ]
df2 = merge(df1, metabolites, by.x=c('eid_1'), by.y=c('Eid')) %>%
  select('eid_1','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL','diag_icd10','icd_group')

Nash = merge(NASH, metabolites, by.x=c('eid'), by.y=c('Eid')) %>%
  select('eid','icd_code','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL')

K760 <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K760'))
K760 <- as.data.frame(K760[!duplicated(K760$eid_1), ])
K760$K760 <- 1
K760 <- select(K760, c(eid_1, K760))

K758 <- subset(hesin_diag, startsWith(as.character(diag_icd10), 'K758'))
K758 <- as.data.frame(K758[!duplicated(K758$eid), ])
K758$K758 <- 1
K758 <- select(K758, c(eid_1, K758))

common5 <- intersect(df2$eid_1, K758$eid_1)
df5_common = df2[common5, ]
K758_common = K758[common5, ]
df5 = merge(df2, K758, by.x=c('eid_1'), by.y=c('eid_1')) %>%
  select('eid_1','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL','diag_icd10') 

MAFLD_df <- subset(df2, startsWith(as.character(diag_icd10), 'K760'))
Fibrosis_cirrhosis <- df2[df2$icd_group == 'K74.0' & !is.na(df2$icd_group), ]
neoplasm <- df2[df2$icd_group == 'C22.0' & !is.na(df2$icd_group), ]
no_liver_disease <- df2[!(df2$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0')), ]

# Nash (K75.8)

NASH_dataframe <- intersect(NASH$eid, metabolites$Eid)
nash_common = NASH[NASH_dataframe, ]
metabolites_common = metabolites[NASH_dataframe, ]
NASH_dataframe = merge(NASH,metabolites, by.x=c('eid'), by.y=c('Eid'))%>% 
  select('eid','L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL')
NASH_df = NASH_dataframe[!duplicated(NASH_dataframe), ]

common4 <- intersect(df2$eid_1, No_Medication$eid)
df4_common = df2[common4, ]
No_Medication_common = No_Medication[common3, ]
healthy_df = merge(df2, No_Medication, by.x=c('eid_1'), by.y=c('eid')) %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL')

healthy_df = healthy_df[!duplicated(healthy_df), ]
metformin_df = metformin_df[!duplicated(metformin_df), ]
count(metformin_df)

## 1.2 Metabolites

common3 <- intersect(df2$eid_1, medication$eid)
df3_common = df2[common3, ]
medication_common = medication[common3, ]
df3 = merge(medication, df2, by.x=c('eid'), by.y=c('eid_1')) %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL','diag_icd10','icd_group')

common4 <- intersect(df3$Eid, covariates$Eid_1)
df3_common = df3[common4, ]
covariates_common = covariates[common4, ]
df4 = merge(covariates, df3, by.x=c('Eid_1'), by.y=c('Eid'))

## 1.3 MRI data 

common5 <- intersect(df4$Eid_1 , mri$eid)
df4_common = df4[common5, ]
mri_common = mri[common5,]
df5 = merge(df3, mri, by.x=c('Eid_1'), by.y=c('eid'))

## 1.5 Death cause

common6 <- intersect(death$eid , death_cause$eid)
death_common = death[common6, ]
death_cause_common = death_cause[common6,]
df6 = merge(death, death_cause, by.x=c('eid_1','ins_index'), by.y=c('eid','ins_index'))

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

# Medications frequent occurrence in Medications  

occurance <- data.frame(table(medication$Medications_2))
occurance2 = occurance[occurance$Freq > 1000,]

# Use mutate_all and recode to transform all row values

for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140875408] <- "allopurinol"} 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140883504] <- "cetirizine" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140916682] <- "evening primrose oil"} 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140910814] <- "sodium thyroxine" } 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140863152] <- "diazepam" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140865580] <- "asacol-400mg-e/c-tablet" } 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871310] <- "ibuprofen" } 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140868226] <- "Aspirin" } 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879616] <- "amitriptyline" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140922174] <- "alendronate_sodium"} 
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140870390] <- "ferrous sulphate" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871310] <- "ibuprofen" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140888366] <- "thiamine_preparation" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140860696] <- "lisinopril" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140864752] <- "lansoprazole" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140871024] <- "vitamin_b_compound_tablet" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879760] <- "bisoprolol" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140879802] <- "amlodipine" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140923346] <- "co-codamol" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141182628] <- "tiotropium" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141176832] <- "seretide_50_evohaler" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1141145660] <- "valsartan" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140863144] <- "zopiclone" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140926606] <- "salbutamol_100micrograms_spacehaler" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140888510] <- "verapamil" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140883066] <- "insulin product" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140860806] <- "ramipril" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140909708] <- "furosemide" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 99999] <- "uncoded" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140864992] <- "tramadol" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] %in% 1140884600] <- "metformin" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140865716] <- "senna" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 2038460150] <- "paracetamol" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140874930] <- "prednisolone" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140888266] <- "warfarin" } 
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140910766] <- "nicorandil" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140874420] <- "quinine" } 
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140879406] <- "ranitidine" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140881856] <- "salbutamol" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141191044] <- "levothyroxine sodium" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 6] <- "Non-cancer Illness"}
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 3004] <- "Mobility problem severity" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140851028] <- "chalk" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140884412] <- "sumatriptan" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140865634] <- "omeprazole" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 2038460150] <- "paracetamol" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140866738] <- "atenolol" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141146234] <- "NA" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140861958] <- "NA" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141194794] <- "bendroflumethiazide" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1141188442] <- "glucosamine product" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140909674] <- "cod_liver_oil_capsule" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140868226] <- "Aspirin" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140852976] <- "multivitamins" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140861998] <- "ventolin 100micrograms inhaler" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1193] <- "omega-3/fish oil supplement" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140909726] <- "vitamin c product" }
for (i in seq_along(medication))  { medication[[i]][medication[[i]] %in% 1140888538] <- "zinc product" }

# Ethnicity values change 

for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1001] <- "British" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 6] <- "Other ethnic group" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3] <- "Asian" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1002] <- "Irish" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3001] <- "Indian" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% -3] <- "Prefer not to answer" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 4001] <- "Caribbean" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1003] <- "Any other white background" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 1002] <- "Irish" }
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 2001] <- "White and Black Caribbean" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3002] <- "Pakistani" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 4002] <- "African" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% 3004] <- "Any other Asian background" } 
for (i in seq_along(metabolites)) { metabolites[[i]][metabolites[[i]] %in% -1] <- "Do not know" }

# Genes changes
for (i in seq_along(genes)) { genes[[i]][genes[[i]] %in% 0] <- "Physical activity" } %>% { genes[[i]][genes[[i]] %in% 1] <- "heterozygous" } %>% { genes[[i]][genes[[i]] %in% 2] <- "homozygous" }

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
hesin_diag[hesin_diag$diag_icd10 %in% c("K75.8"), "icd_group"] <- "NASH"

# ICD CODE grouping

death_cause[death_cause$cause_icd10 %in% c("K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "icd_group"] <- "K70.0"
death_cause[death_cause$cause_icd10 %in% c("K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718","K719"), "icd_group"] <- "K71.0"
death_cause[death_cause$cause_icd10 %in% c("K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728","K729"), "icd_group"] <- "K72.0"
death_cause[death_cause$cause_icd10 %in% c("K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "icd_group"] <- "K73.0"
death_cause[death_cause$cause_icd10 %in% c("K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "icd_group"] <- "K74.0"
death_cause[death_cause$cause_icd10 %in% c("K751", "K752", "K753", "K754", "K755", "K756", "K757"), "icd_group"] <- "K75.0"
death_cause[death_cause$cause_icd10 %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767"), "icd_group"] <- "K76.0"
death_cause[death_cause$cause_icd10 %in% c("K771", "K772", "K773", "K774", "K705", "K776", "K777"), "icd_group"] <- "K77.0"
death_cause[death_cause$cause_icd10 %in% c("C221", "C222", "C223", "C224", "C225", "C226", "C227"), "icd_group"] <- "c22.0"
death_cause[death_cause$cause_icd10 %in% c("K75.8"), "icd_group"] <- "NASH"

death_cause[death_cause$icd_group %in% c("K70.0"), "icd_names"] <- "Alcoholic liver disease"
death_cause[death_cause$icd_group %in% c("K71.0"), "icd_names"] <- "Toxic liver disease "
death_cause[death_cause$icd_group %in% c("K72.0","K73.0"), "icd_names"] <- "Hepatic failure"
death_cause[death_cause$icd_group %in% c("K74.0"), "icd_names"] <- "Fibrosis and cirrhosis"
death_cause[death_cause$icd_group %in% c("K75.0"), "icd_names"] <- "Inflammatory liver diseases"
death_cause[death_cause$icd_group %in% c("K75.8"), "icd_names"] <- "NASH"
death_cause[death_cause$icd_group %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767"), "icd_names"] <- "Chronic Hepatitis MAFLD"
death_cause[death_cause$icd_group %in% c("C22.0"), "icd_names"] <- "liver Malignant neoplasm"


baseline_df[baseline_df$Diagnosis_grp %in% c("K700", "K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "Diagnosis_grp"] <- "K70.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K710", "K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718", "K719"), "Diagnosis_grp"] <- "K71.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K720", "K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728", "K729"), "Diagnosis_grp"] <- "K72.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K730", "K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "Diagnosis_grp"] <- "K73.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K740", "K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "Diagnosis_grp"] <- "K74.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K750", "K751", "K752", "K753", "K754", "K755", "K756", "K757", "K759"), "Diagnosis_grp"] <- "K75.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K760", "K761", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769"), "Diagnosis_grp"] <- "K76.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("K770", "K771", "K772", "K773", "K774", "K705", "K776", "K777", "K778", "K779"), "Diagnosis_grp"] <- "K77.0"
baseline_df[baseline_df$Diagnosis_grp %in% c("C220", "C221", "C222", "C223", "C224", "C705", "C226", "C227", "C228", "C229"), "Diagnosis_grp"] <- "C22.0"