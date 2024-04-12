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
library(plotly)
library(readr)
library(ggrepel)
library(EnhancedVolcano)
library(effectsize)
library(broom)
library(MatchIt)
library(cobalt)

# Read Text and CSV Files_the R Script

medication <- read.csv(file = 'C:/Users/User/Desktop/Data/medication.csv')[ ,2:21]
covariates2 <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates_Markus.csv') %>%
  select('eid','X52.0.0')
hesin_diag <- read.table(file = 'C:/Users/User/Desktop/Data/hesin_diag.txt' ,header = TRUE , sep = "\t" )%>% 
  select('requires.eid','ins_index','diag_icd10')
hesin <- read.table(file = 'C:/Users/User/Desktop/Data/hesin.txt', header = TRUE , sep = "\t") %>% 
  select('eid','ins_index','epistart','epiend','epidur','disdate','admidate')
metabolites <- read.csv(file = 'C:/Users/User/Desktop/Data/metabolites.csv')
mri <- read.csv(file = 'C:/Users/User/Desktop/Data/MRI.csv')%>%
  select('eid','X40061.2.0')
death <- read.table(file = 'C:/Users/User/Desktop/Data/Endpoints/death/death.txt', header = TRUE , sep = "\t")%>% 
  select('eid','ins_index','date_of_death')
death_cause <- read.table(file = 'C:/Users/User/Desktop/Data/Endpoints/death/death_cause.txt', header = TRUE , sep = "\t")%>% 
  select('eid','ins_index','cause_icd10')
covariates <- read.csv(file = 'C:/Users/User/Desktop/Data/covariates.csv')
genes <- read.csv(file = 'C:/Users/User/Desktop/Data/genes.csv')  %>%
  select("eid","rs2642438_A","rs1260326_T","rs780094_T","rs72613567_TA","rs10069690_T","rs9272105_A","rs2856718_T","rs7453920_A","rs3077_G",
         "rs9277535_G","rs58542926_T","rs455804_A","rs738409_G")  
NASH <- read.csv(file = 'C:/Users/User/Desktop/Data/Results/Nash.csv')

# After setnames of metabolites:

metabolite <- metabolites %>% select('Eid','Ethnicity','BMI','Age','Sex')

# 'Total_Free_C','Clinical_LDL_C','HDL_C','Total_Triglycerides','Linoleic_Acid','Omega-3 _FA_','Omega-6__FA_','Triglycerides_VLDL', 'Triglycerides_LDL','Triglycerides_HDL','CE_CM and EL_VLDL','CE_VL_VLDL','CE_L_VLDL','CE_M_VLDL','CE_S_VLDL','CE_VS_VLDL','	CE_IDL','CE_L_LDL','CE_M_LDL','CE_S_LDL','CE_VL_HDL','CE_L_HDL','CE_L_HDL',

# Mutation of data for Medication
 
setnames(hesin_diag,"requires.eid",'eid')
setnames(medication, old=colnames(medication), new = c('eid','Medications_1','Medications_2','Medications_3','Medications_4','Medications_5','Medications_6','Medications_7','Medications_8','Medications_9','Medications_10','Medications_11','Medications_12','Medications_13','Medications_14','Medications_15','Medications_16','Medications_17','Medications_18','Medications_19'))
setnames(metabolites, old=colnames(metabolites), new = c('No.','eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
                                                         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
                                                         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
                                                         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
                                                         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','Ethnicity','BMI','Age','Sex'))
setnames(mri,old=colnames(mri), new = c('eid',"Liver_fat"))
setnames(covariates,old=colnames(covariates),new = c('eid','ALT','AST','GGT'))
setnames(covariates2,old=colnames(covariates2),new = c('eid','init_exam'))
setnames(death,"eid","eid_1")
setnames(NASH,"X41202.0.0","icd_code")
setnames(K760,"K760","Diagnosis")

#'TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLC','PC','S','AB','A1','ABA1','TLFA','DU','O3FA'

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

admissionP <- na.omit(hesin) %>% select("eid","admidate")

# death info consists several columns mentioned below for MAFLD and NASH associated death.  

death_info <- na.omit(df6) %>% select("eid_1","date_of_death", "icd_group","icd_names")

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
         
# No medication Dataframe generated_calculations.R file with no usage of medications for NASH and MAFLD patients. 

common4 <- intersect(df2$eid_1, No_Medication$eid)
df4_common = df2[common4, ]
No_Medication_common = No_Medication[common3, ]
healthy_df = merge(df2, No_Medication, by.x=c('eid_1'), by.y=c('eid')) %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL')

healthy_df = healthy_df[!duplicated(healthy_df), ]
healthy_metformin_df <- metformin_df[!(metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0')), ]
diagnosed_metformin_df <- metformin_df[metformin_df$icd_group %in% c('K70.0', 'K71.0', 'K72.0', 'K73.0', 'K74.0', 'K75.0', 'K76.0', 'K77.0', 'C22.0'), ]
MAFLD_df <- diagnosed_metformin_df[diagnosed_metformin_df$icd_group == 'K76.0', ]


## 1.5 Death cause and death dataframe intersect to observe cause of NAFLD and MAFLD. 

common6 <- intersect(death$eid_1 , death_cause$eid)
death_common = death[common6, ]
death_cause_common = death_cause[common6,]
death_df = merge(death, death_cause, by.x=c('eid_1','ins_index'), by.y=c('eid','ins_index'))
death_df = death_df[!duplicated(death_df$eid_1), ]

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
for (i in seq_along(medication)) { medication[[i]][medication[[i]] ==1140860806] <- "ramipril" }
for (i in seq_along(medication)) { medication[[i]][medication[[i]] == 1140884600] <- "metformin"}
for (i in seq_along(medication)) { medication[[i]][medication[[i]] == 1141171646] <- "pioglitazone"}

# Ethnicity values change 

# Genes changes
for (i in seq_along(genes)) { genes[[i]][genes[[i]] %in% 0] <- "Physical activity" } %>% { genes[[i]][genes[[i]] %in% 1] <- "heterozygous" } %>% { genes[[i]][genes[[i]] %in% 2] <- "homozygous" }

# Hesin diagnosis coding :

hesin_diag[hesin_diag$diag_icd10 %in% c("K700", "K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "icd_group"] <- "K70.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K710", "K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718", "K719"), "icd_group"] <- "K71.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K720", "K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728", "K729"), "icd_group"] <- "K72.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K730", "K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "icd_group"] <- "K73.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K740", "K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "icd_group"] <- "K74.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K750", "K751", "K752", "K753", "K754", "K755", "K756", "K757", "K759"), "icd_group"] <- "K75.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769"), "icd_group"] <- "K76.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K770", "K771", "K772", "K773", "K774", "K705", "K776", "K777", "K778", "K779"), "icd_group"] <- "K77.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K220", "K221", "K222", "K223", "K224", "K705", "K226", "K227", "K228", "K229"), "icd_group"] <- "C22.0"
hesin_diag[hesin_diag$diag_icd10 %in% c("K75.8"), "icd_group"] <- "K75.8"
hesin_diag[hesin_diag$diag_icd10 %in% c("K76.0"), "icd_group"] <- "K76.0"

# ICD CODE grouping

death_cause[death_cause$cause_icd10 %in% c("K701", "K702", "K703", "K704", "K705", "K706", "K707", "K708", "K709"), "icd_group"] <- "K70.0"
death_cause[death_cause$cause_icd10 %in% c("K711", "K712", "K713", "K714", "K715", "K716", "K717", "K718","K719"), "icd_group"] <- "K71.0"
death_cause[death_cause$cause_icd10 %in% c("K721", "K722", "K723", "K724", "K725", "K726", "K727", "K728","K729"), "icd_group"] <- "K72.0"
death_cause[death_cause$cause_icd10 %in% c("K731", "K732", "K733", "K734", "K735", "K736", "K737", "K738", "K739"), "icd_group"] <- "K73.0"
death_cause[death_cause$cause_icd10 %in% c("K741", "K742", "K743", "K744", "K745", "K746", "K747", "K748", "K749"), "icd_group"] <- "K74.0"
death_cause[death_cause$cause_icd10 %in% c("K751", "K752", "K753", "K754", "K755", "K756", "K757", "K758", "K759"), "icd_group"] <- "K75.0"
death_cause[death_cause$cause_icd10 %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767","K768", "K769"), "icd_group"] <- "K76.0"
death_cause[death_cause$cause_icd10 %in% c("K771", "K772", "K773", "K774", "K705", "K776", "K777","K778", "K779"), "icd_group"] <- "K77.0"
death_cause[death_cause$cause_icd10 %in% c("C221", "C222", "C223", "C224", "C225", "C226", "C227","C228", "C229"), "icd_group"] <- "C22.0"

death_cause[death_cause$icd_group %in% c("K70.0"), "icd_names"] <- "Alcoholic liver disease"
death_cause[death_cause$icd_group %in% c("K71.0"), "icd_names"] <- "Toxic liver disease "
death_cause[death_cause$icd_group %in% c("K72.0","K73.0"), "icd_names"] <- "Hepatic failure"
death_cause[death_cause$icd_group %in% c("K74.0"), "icd_names"] <- "Fibrosis and cirrhosis"
death_cause[death_cause$icd_group %in% c("K75.0"), "icd_names"] <- "Inflammatory liver diseases"
death_cause[death_cause$icd_group %in% c("K75.8"), "icd_names"] <- "NASH"
death_cause[death_cause$icd_group %in% c("K761", "K762", "K763", "K764", "K765", "K766", "K767"), "icd_names"] <- "Chronic Hepatitis MAFLD"
death_cause[death_cause$icd_group %in% c("C22.0"), "icd_names"] <- "liver Malignant neoplasm"
