# Libraries to support different functions in this file. 

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(lubridate)
library(plyr)
library(readr)
library(janitor)
library(DT)
library(vtree)
library(Matrix)
library(ff)

df22 <- read.csv(file = 'C:/Users/User/Desktop/Data/metabolites.csv')
data3 <- read.csv(file = 'C:/Users/User/Desktop/Data/MRI.csv')
setnames(df22, old=colnames(df22), new = c('No.','Eid','Total_Cholesterol','Total_Cholesterol-HDL_C','Remnant_Cholesterol','VLDL_Cholesterol','Clinical_LDL_Cholesterol','LDL_Cholesterol','HDL_Cholesterol','Total_Triglycerides','Triglycerides_in_VLDL', 'Triglycerides_in_LDL','Triglycerides_in_HDL','Total_Phospholipids_in_Lipoprotein_Particles','Phospholipids_in_VLDL', 'Phospholipids_in_LDL','Phospholipids_in_HDL','Total_Esterified_Cholesterol', 'Cholesteryl_Esters_in_VLDL','Cholesteryl_Esters_in_LDL','Cholesteryl_Esters_in_HDL','Total_Free_Cholesterol','Free_Cholesterol_in_VLDL','Free_Cholesterol_in_LDL','Free Cholesterol in HDL', 'Total Lipids in Lipoprotein Particles', 'Total Lipids in VLDL','Total Lipids in LDL','Total Lipids in HDL','Total Concentration of Lipoprotein Particles','Concentration of VLDL Particles','Concentration of LDL Particles',	'Concentration_of_HDL_Particles','Average Diameter for VLDL Particles','Average Diameter for LDL Particles,Average Diameter for HDL_Particles','Phosphoglycerides','Triglycerides_to_Phosphoglycerides ratio','Total_Cholines','Phosphatidylcholines','Sphingomyelins','Apolipoprotein_B','Apolipoprotein_A1','Apolipo protein B-Apolipoprotein A1 ratio',
                                           'Total_Fatty_Acids','Degree of Unsaturation','Omega-3 Fatty_Acids','Omega-6_Fatty_Acids','Polyunsaturated_Fatty_Acids','Monounsaturated_Fatty_Acids','Saturated_Fatty_Acids','Linoleic_Acid','Docosahexaenoic_Acid','Omega-3_Fatty_Acids_to_Total_Fatty_Acids_percentage','Omega-6 Fatty Acids to Total Fatty Acids percentage','Polyunsaturated Fatty Acids to Total Fatty Acids percentage','Monounsaturated Fatty Acids to Total Fatty Acids percentage','Saturated Fatty Acids to Total Fatty Acids percentage','Linoleic Acid to Total Fatty Acids percentage','Docosahexaenoic Acid to Total Fatty Acids percentage','Polyunsaturated Fatty Acids to Monounsaturated Fatty Acids ratio','Omega-6 Fatty Acids to Omega-3 Fatty Acids ratio','Alanine','Glutamine','Glycine','Histidine','Total Concentration of Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	Concentration of Chylomicrons and EL_VLDL Particles','Total Lipids in Chylomicrons and EL_VLDL',
                                           'Phospholipids in Chylomicrons and EL_VLDL','Cholesterol in Chylomicrons and EL_VLDL','Cholesteryl Esters in Chylomicrons and EL_VLDL','Free Cholesterol in Chylomicrons and EL_VLDL','Triglycerides in Chylomicrons and EL_VLDL','Concentration of VL_VLDL Particle','Total Lipids in VL_VLDL','Phospholipids in VL_VLD','Cholesterol in VL_VLDL','Cholesteryl Esters in VL_VLDL','	Free Cholesterol in VL_VLDL','Triglycerides in VL_VLDL','Concentration of L_VLDL Particles','Total Lipids in L_VLDL','Phospholipids in L_VLDL','Cholesterol in L_VLDL','Cholesteryl Esters in L_VLDL','Free Cholesterol in L_VLDL','Triglycerides in L_VLDL','Concentration of M_VLDL','Total Lipids in M_VLDL','Phospholipids in M_VLDL','Cholesterol in M_VLDL','Cholesteryl Esters in M_VLDL','Free Cholesterol in M_VLDL','Triglycerides in M_VLDL','	Concentration of S_VLDL Particles','Total Lipids in S_VLDL','Phospholipids in S_VLDL','Cholesterol in S_VLDL','Cholesteryl Esters in S_VLDL','Free Cholesterol in S_VLDL','Triglycerides in S_VLDL',
                                           'Concentration of VS_VLDL Particles','Total Lipids in VS_VLDL','Phospholipids in_VS_VLDL','Cholesterol in_VS_VLDL','Cholesteryl Esters in VS_VLDL','Free Cholesterol in VS_VLDL','Triglycerides in VS_VLDL','Concentration of IDL Particles','Total Lipids in IDL','Phospholipids in IDL','Cholesterol in IDL','	Cholesteryl Esters in IDL','Free_Cholesterol_in_IDLs','Triglycerides in IDL','Concentration of L_LDL Particles','	Total Lipids in L_LDL','Phospholipids in L_LDL','Cholesterol in L_LDL','Cholesteryl Esters in L_LDL','Free Cholesterol in L_LDL','Triglycerides in L_LDL','Concentration of M_LDL Particles','Total Lipids in M_LDL','Phospholipids in M_LDL','Cholesterol in M_LDL','Cholesteryl Esters in M_LDL','Free Cholesterol in M_LDL','Triglycerides in M_LDL','Concentration of S_LDL Particles','Total Lipids in S_LDL','Phospholipids in S_LDL','Cholesterol in S_LDL','Cholesteryl Esters in S_LDL','Free Cholesterol in S_LDL','Triglycerides in S_LDL','Concentration of VL_HDL_Particles','Total Lipids in VL_HDL','Phospholipids in VL_HDL','Cholesterol in VL_HDL',
                                           'Cholesteryl Esters in VL_HDL','Free Cholesterol in VL_HDL','Triglycerides in_VL_HDL','Concentration of L_HDL_Particles','Total Lipids in L_HDL','Phospholipids in L_HDL','Cholesterol in L_HDL','Cholesteryl Esters in L_HDL','Free Cholesterol in L_HDL','Triglycerides in L_HDL','Concentration of M_HDL_Particles','Total Lipids in M_HDL','Phospholipids in M_HDL','Cholesterol in M_HDL','Ethnicity','BMI','Age','Sex'))
setnames(data3,old=colnames(data3), new = c('X','Id',"Liver_fat"))

## 1.4 Metabolites 

common4 <- intersect(df4$Eid, df22$Eid)
df4_common = df4[common3, ]
df22_common = df22[common3, ]
common4_comb = merge(df4, df22, by.x=c('Eid'), by.y=c('Eid'))
df5 = common4_comb

## 1.5 df6 consists with Patients columns (1-19), with basic diagnosis code   

common5 <- intersect(df5$Eid, data3$Id)
df5_common = df5[common4, ]
data3_common = data3[common4, ]
common5_comb = merge(df5, data3, by.x=c('Eid'), by.y=c('Id'))%>% 
  select('Eid','Patient_1', 'Patient_2','Patient_3','Patient_4','Patient_5','Patient_6','Patient_7','Patient_8','Patient_9','Patient_10','Patient_11','Patient_12','Patient_13','Patient_14','Patient_15','Patient_16','Patient_17','Patient_18','Diag-code','epistart','epiend', 'admidate', 'disdate','Gender','Birth_Year','Liver_fat','Total_Cholesterol','Total_Free_Cholesterol','Total_Triglycerides','LDL_Cholesterol','HDL_Cholesterol','Total_Cholesterol-HDL_C')
df6 = common5_comb

