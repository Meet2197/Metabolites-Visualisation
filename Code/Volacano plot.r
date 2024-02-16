# Library to upload : 

library(PheWAS)
library(broom)
library(ggplot2)
library(ggrepel)
library(effectsize)
library(plotly)


# File upload : 

metabolite_metformin <-merge(ALL, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','metformin')

# summary table loop :

summary_table <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_metformin)
    first_tables[[i]] <- broom::tidy(models[[i]])
    
    
  }  
  
  # Combine the rows together into a final table
  final_table <- do.call("rbind", first_tables)
  
  return(final_table)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table1 <- summary_table(metabolite_metformin[2:170])

# metformin filter with table : 

final_table_Metformin <- subset(final_table1, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_table_Metformin$logp <- -(log10(final_table_Metformin$p.value))
final_table_Metformin <- final_table_Metformin[-c(1,2,3,4,5)]  

# generate smd value. 

smd_values <- apply(metabolite_metformin[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_metformin$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_metformin$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_df <- data.frame(smd_values)

# Adds new row to df

smd_values[3,]=c(3,smd_values[2,-1]/smd_values[1,-1])

# Add another new row at index 4
smd_df[4,] <- c(4, log2(smd_df[3,]))

# Remove the first three rows
table_smd<-as.data.frame(t(as.data.frame(smd_df)))
table_smd<-table_smd[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_Metformin<-cbind(table_smd, final_table_Metformin)

# name change for smd 
colnames(all_Metformin)[colnames(all_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_Metformin)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of Metformin Data", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# Generate volcano plot

Diabetes_metformin <- ALL %>%
  mutate(Diabetes = ifelse(Diabetes %in% c(-3, -1, 0), 0, Diabetes)) %>%
  filter(!is.na(Diabetes))

# Create a new variable to count occurrences
Diabetes_metformin1 <- Diabetes_metformin %>%
  group_by(Diabetes, metformin) %>%
  summarise(Count = n())

# Convert factors to meaningful labels
Diabetes_metformin1$Diabetes <- ifelse(Diabetes_metformin1$Diabetes == 0, "Patients without Diabetes", "Patients with Diabetes")
Diabetes_metformin1$metformin <- ifelse(Diabetes_metformin1$metformin == 0, "0", "1")

# Generate the plot
ggplot(Diabetes_metformin1, aes(x = Diabetes, y = Count, fill = metformin)) +
  geom_col(position = "dodge") + # Use geom_col() instead of geom_bar()
  geom_text(aes(label = Count), position = position_dodge(width = 1), vjust = -0.5, color = "black") +
  labs(title = "Comparison of Metformin Users vs Non-Users among Diabetes Patients",
       x = "Diabetes Status",
       y = "Count") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "red"), labels = c("Non-Users", "Users")) + # Specify fill colors for metformin
  theme_minimal()


# Diabetic Patients 

metabolite_Diabetic <-merge(ALL, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','metformin','Diabetes')

metabolite_Diabetic <- metabolite_Diabetic %>%
  filter(Diabetes == 1)

summary <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_Diabetic)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table <- do.call("rbind", first_tables)
  
  return(final_table)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table <- summary(metabolite_Diabetic[2:170])

# metformin filter with table : 

final_table_Diabetes <- subset(final_table, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_table_Diabetes$logp <- -(log10(final_table_Diabetes$p.value))
final_table_Diabetes <- final_table_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_diabetes <- apply(metabolite_Diabetic[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_Diabetic$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_Diabetic$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_Diabetic_df <- data.frame(smd_diabetes)

# Adds new row to df

smd_Diabetic_df[3,]=c(3,smd_Diabetic_df[2,-1]/smd_Diabetic_df[1,-1])

# Add another new row at index 4
smd_Diabetic_df[4,] <- c(4, log2(smd_Diabetic_df[3,]))

# Remove the first three rows
table_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_Diabetic_df)))
table_smd_Diabetic<-table_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_diabetes_Metformin<-cbind(table_smd_Diabetic, final_table_Diabetes)

# name change for smd 
colnames(all_diabetes_Metformin)[colnames(all_diabetes_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_diabetes_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_diabetes_Metformin)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of Metformin Data in Diabitic", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# nash Patients 

metabolite_nash <-merge(ALL, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','metformin','nash')

metabolite_nash <- metabolite_nash %>%
  filter(nash == 1)

summary2 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_nash)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table <- do.call("rbind", first_tables)
  
  return(final_table)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table2 <- summary2(metabolite_nash[2:170])

# metformin filter with table : 

final_table_nash <- subset(final_table2, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_table_nash$logp <- -(log10(final_table_nash$p.value))
final_table_nash <- final_table_nash[-c(1,2,3,4,5)]  

# generate smd value. 

smd_nash <- apply(metabolite_nash[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_nash$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_nash$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_nash_df <- data.frame(smd_nash)

# Adds new row to df

smd_nash_df[3,]=c(3,smd_nash_df[2,-1]/smd_nash_df[1,-1])

# Add another new row at index 4
smd_nash_df[4,] <- c(4, log2(smd_nash_df[3,]))

# Remove the first three rows
t_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_nash_df)))
t_smd_Diabetic<-t_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_nash_Metformin<-cbind(t_smd_Diabetic, final_table_nash)

# name change for smd 
colnames(all_nash_Metformin)[colnames(all_nash_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_nash_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_nash_Metformin)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of nash Patients consuming Metformin", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

metabolite_nash_Diabetes <-merge(ALL, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','metformin','nash','Diabetes')

metabolite_nash_Diabetes <- metabolite_nash_Diabetes %>%
  filter(nash == 1, Diabetes == 1)

summary3 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_nash_Diabetes)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table3 <- do.call("rbind", first_tables)
  
  return(final_table3)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table3 <- summary3(metabolite_nash_Diabetes[2:170])

# metformin filter with table : 

final_nash_Diabetes <- subset(final_table3, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_nash_Diabetes$logp <- -(log10(final_nash_Diabetes$p.value))
final_nash_Diabetes <- final_nash_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_nash_diabetes <- apply(metabolite_nash_Diabetes[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_nash_Diabetes$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_nash_Diabetes$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_nash_diabetes_df <- data.frame(smd_nash_diabetes)

# Adds new row to df

smd_nash_diabetes_df[3,]=c(3,smd_nash_diabetes_df[2,-1]/smd_nash_diabetes_df[1,-1])

# Add another new row at index 4
smd_nash_diabetes_df[4,] <- c(4, log2(smd_nash_diabetes_df[3,]))

# Remove the first three rows
t_smd_nash_Diabetic<-as.data.frame(t(as.data.frame(smd_nash_diabetes_df)))
t_smd_nash_Diabetic<-t_smd_nash_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_nash_Metformin_Diabetic<-cbind(t_smd_nash_Diabetic, final_nash_Diabetes)

# name change for smd 
colnames(all_nash_Metformin)[colnames(all_nash_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_nash_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_nash_Metformin)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of nash Patients consuming Metformin(Diabetic)", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# pioglitazone

metabolite_pioglitazone <-merge(ALL_pioglitazone, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','pioglitazone')

# summary table loop for pioglitazone:

summary_table4 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = metabolite_pioglitazone)
    first_tables[[i]] <- broom::tidy(models[[i]])
    
    
  }  
  
  # Combine the rows together into a final table
  final_table4 <- do.call("rbind", first_tables)
  
  return(final_table4)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 

final_table4 <- summary_table4(metabolite_pioglitazone[2:170])

# pioglitazone filter with table : 

final_table_pioglitazone <- subset(final_table4, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

final_table_pioglitazone$logp <- -(log10(final_table_pioglitazone$p.value))
final_table_pioglitazone <- final_table_pioglitazone[-c(1,2,3,4,5)]  

# generate smd value. 

smd_values_pioglitazone <- apply(metabolite_pioglitazone[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_pioglitazone$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_pioglitazone$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_pioglitazone <- data.frame(smd_values_pioglitazone)

# Adds new row to df

smd_pioglitazone[3,]=c(3,smd_pioglitazone[2,-1]/smd_pioglitazone[1,-1])

# Add another new row at index 4
smd_pioglitazone[4,] <- c(4, log2(smd_pioglitazone[3,]))

# Remove the first three rows
table_smd_pioglitazone<-as.data.frame(t(as.data.frame(smd_pioglitazone)))
table_smd_pioglitazone<-table_smd_pioglitazone[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
all_pioglitazone<-cbind(table_smd_pioglitazone, final_table_pioglitazone)

# name change for smd 
colnames(all_pioglitazone)[colnames(all_pioglitazone) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_pioglitazone, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_pioglitazone)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of pioglitazone Data", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# Diabetic Patients 

Diabetic_pioglitazone <-merge(ALL_pioglitazone, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','pioglitazone','Diabetes')

Diabetic_pioglitazone <- Diabetic_pioglitazone %>%
  filter(Diabetes == 1)

summary5 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = Diabetic_pioglitazone)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table5 <- do.call("rbind", first_tables)
  
  return(final_table5)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 

final_table5 <- summary5(Diabetic_pioglitazone[2:170])

# pioglitazone filter with table : 

pioglitazone_Diabetes <- subset(final_table5, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

pioglitazone_Diabetes$logp <- -(log10(pioglitazone_Diabetes$p.value))
pioglitazone_Diabetes <- pioglitazone_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_diabetes <- apply(Diabetic_pioglitazone[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(Diabetic_pioglitazone$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(Diabetic_pioglitazone$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_Diabetic_df <- data.frame(smd_diabetes)

# Adds new row to df

smd_Diabetic_df[3,]=c(3,smd_Diabetic_df[2,-1]/smd_Diabetic_df[1,-1])

# Add another new row at index 4
smd_Diabetic_df[4,] <- c(4, log2(smd_Diabetic_df[3,]))

# Remove the first three rows
table_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_Diabetic_df)))
table_smd_Diabetic<-table_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
all_diabetes_pioglitazone<-cbind(table_smd_Diabetic, pioglitazone_Diabetes)

# name change for smd 
colnames(all_diabetes_pioglitazone)[colnames(all_diabetes_pioglitazone) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_diabetes_pioglitazone, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_diabetes_pioglitazone)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of pioglitazone for Diabitic Patient", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# MASLD Patients with pioglitazone :

pioglitazone_MASLD <-merge(ALL_pioglitazone, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','pioglitazone','MASLD')

pioglitazone_MASLD <- pioglitazone_MASLD %>%
  filter(MASLD == 1)

summary8 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = pioglitazone_MASLD)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table8 <- do.call("rbind", first_tables)
  
  return(final_table8)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 

final_table8 <- summary8(pioglitazone_MASLD[2:170])

# pioglitazone filter with table : 

final_table_MASLD <- subset(final_table8, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

final_table_MASLD$logp <- -(log10(final_table_MASLD$p.value))
final_table_MASLD1 <- final_table_MASLD[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD_pioglitazone <- apply(pioglitazone_MASLD[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(pioglitazone_MASLD$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(pioglitazone_MASLD$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASLD_pioglitazone <- data.frame(smd_MASLD_pioglitazone)

# Adds new row to df

smd_MASLD_pioglitazone[3,]=c(3,smd_MASLD_pioglitazone[2,-1]/smd_MASLD_pioglitazone[1,-1])

# Add another new row at index 4
smd_MASLD_pioglitazone[4,] <- c(4, log2(smd_MASLD_pioglitazone[3,]))

# Remove the first three rows
t_smd_MASLD_pioglitazone<-as.data.frame(t(as.data.frame(smd_MASLD_pioglitazone)))
t_smd_MASLD_pioglitazone<-t_smd_MASLD_pioglitazone[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
t_smd_MASLD_pioglitazone<-cbind(t_smd_MASLD_pioglitazone, final_table_MASLD1)

# name change for smd 
colnames(t_smd_MASLD_pioglitazone)[colnames(t_smd_MASLD_pioglitazone) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 
ggplot(t_smd_MASLD_pioglitazone, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(t_smd_MASLD_pioglitazone)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of MASLD Patients consuming pioglitazone", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

# MASLD diabetic Patients consuming pioglitazone :

metabolite_MASLD_Diabetes <-merge(ALL_pioglitazone, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','Tl_C','Tl_C-HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Tl_TG','TG_VLDL', 'TG_LDL','TG_HDL','Tl_P_Lipoprotein','P_VLDL', 'P_LDL','P_HDL','Tl_Esterified_C', 'CE_VLDL','CE_LDL','CE_HDL','Tl_Free_C','FC_VLDL','FC_LDL','FC_HDL', 'Tl_LP_Lipoprotein ', 'Tl_LP_VLDL','Tl_LP_LDL','Tl_LP_HDL','Tl  Lipoprotein ','VLDL','LDL','Conc._of_HDL','Avg. Diameter for VLDL ','Avg. Diameter for LDL ,Avg. Diameter for HDL_','Pho','TG_to_Phosphoglycerides ratio','Tl_Cholines','Phosphatidylcholines','Sphingomyelins','Apo_B','Apo_A1','Apolipo protein B-Apo A1 ratio',
         'Tl_FA_','Degree of Unsaturation','Omega-3_FA_','Omega-6_FA_','Polyunsaturated__FA_','Monounsaturated__FA_','Saturated__FA_','Linoleic_Acid','DHA','Omega-3_FA_to_Tl_FA_Per','Omega-6_FA_to Tl_FA_Per','Polyunsaturated_FA_to Tl_FA_Per','Monounsaturated_FA_to Tl_FA_Per','Saturated_FA_to Tl_FA_Per','LA_to Tl_FA_Per','Docosahexaenoic Acid to Tl_FA_Per','Polyunsaturated_FA_to Monounsaturated_FA_ratio','Omega-6_FA_to Omega-3_FA_ratio','Alanine','Glutamine','Glycine','Histidine','Tl  Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)','Isoleucine','Leucine','Valine','Phenylalanine','Tyrosine','Glucose','Lactate','Pyruvate','Citrate','3-Hydroxybutyrate','Acetate','Acetoacetate','Acetone','Creatinine','Albumin','Glycoprotein Acetyls','	 CM and EL_VLDL ','Tl_LP_CM and EL_VLDL',
         'P_CM and EL_VLDL','C_CM and EL_VLDL','CE_CM and EL_VLDL','FC_CM and EL_VLDL','TG_CM and EL_VLDL',' VL_VLDL Particle','Tl_LP_VL_VLDL','P_VL_VLD','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL','Tl_LP_L_VLDL','P_L_VLDL','C_L_VLDL','CE_L_VLDL','FC_L_VLDL','TG_L_VLDL','M_VLDL','Tl_LP_M_VLDL','P_M_VLDL','C_M_VLDL','CE_M_VLDL','FC_M_VLDL','TG_M_VLDL','S_VLDL','Tl_LP_S_VLDL','P_S_VLDL','C_S_VLDL','CE_S_VLDL','FC_S_VLDL','TG_S_VLDL',
         'VS_VLDL','Tl_LP_VS_VLDL','Ph__VS_VLDL','C_VS_VLDL','CE_VS_VLDL','FC_VS_VLDL','TG_VS_VLDL','IDL','Tl_LP_IDL','P_IDL','C_IDL','CE_IDL','FC_IDLs','TG_IDL','L_LDL','Tl_LP_L_LDL','P_L_LDL','C_L_LDL','CE_L_LDL','FC_L_LDL','TG_L_LDL','M_LDL','Tl_LP_M_LDL','P_M_LDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','Tl_LP_S_LDL','P_S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL',' VL_HDL','Tl_LP_VL_HDL','P_VL_HDL','C_VL_HDL',
         'CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','Tl_LP_L_HDL','P_L_HDL','C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','Tl_LP_M_HDL','P_M_HDL','C_M_HDL','pioglitazone','MASLD','Diabetes')

# MASLD and diabetes filter:

metabolite_MASLD_Diabetes <- metabolite_MASLD_Diabetes %>%
  filter(MASLD == 1, Diabetes == 1)

#Summary generation :

summary9 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = metabolite_MASLD_Diabetes)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table9 <- do.call("rbind", first_tables)
  
  return(final_table9)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 

final_table9 <- summary9(metabolite_MASLD_Diabetes[2:170])

# pioglitazone filter with table : 

final_MASLD_Diabetes <- subset(final_table9, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

final_MASLD_Diabetes$logp <- -(log10(final_MASLD_Diabetes$p.value))
final_MASLD_Diabetes <- final_MASLD_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD_diabetes <- apply(metabolite_MASLD_Diabetes[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_MASLD_Diabetes$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_MASLD_Diabetes$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASLD_diabetes_df <- data.frame(smd_MASLD_diabetes)

# Adds new row to df

smd_MASLD_diabetes_df[3,]=c(3,smd_MASLD_diabetes_df[2,-1]/smd_MASLD_diabetes_df[1,-1])

# Add another new row at index 4
smd_MASLD_diabetes_df[4,] <- c(4, log2(smd_MASLD_diabetes_df[3,]))

# Remove the first three rows
t_smd_MASLD_Diabetic<-as.data.frame(t(as.data.frame(smd_MASLD_diabetes_df)))
t_smd_MASLD_Diabetic<-t_smd_MASLD_Diabetic[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
all_MASLD_pioglitazone_Diabetic<-cbind(t_smd_MASLD_Diabetic, final_MASLD_Diabetes)

# name change for smd 
colnames(all_MASLD_pioglitazone_Diabetic)[colnames(all_MASLD_pioglitazone_Diabetic) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 
ggplot(all_MASLD_pioglitazone_Diabetic, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) + # Color points by significance
  geom_text(aes(label = rownames(all_MASLD_pioglitazone_Diabetic)), size = 3, vjust = -0.5, hjust = 0.5) + # Add row names
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of MASLD Patients consuming pioglitazone(Diabetic)", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() #  minimal theme for aesthetics

