
# Dataset from NASH_df, MAFLD and healthy Patients dataframes merged with sub metabolites mentioned.  

combined_data <- rbind(
  mutate(Nash, Data_Frame = "NASH_df"),
  mutate(MAFLD, Data_Frame = "MAFLD"),
  mutate(Healthy_Patient, Data_Frame = "healthy_df"))

# Create the scatter plot with different colors for different variables with ggplot. 

ggplot(combined_data, aes(x = Data_Frame, y = Re_C, color = "Re_C")) +
  geom_jitter(width = 0.2) +
  geom_point(data = combined_data, aes(x = Data_Frame, y = Total_C, color = "Tl_C")) +
  geom_point(data = combined_data, aes(x = Data_Frame, y = VLDL_C, color = "VLDL_C")) +
  geom_point(data = combined_data, aes(x = Data_Frame, y = Cl_LDL_C, color = "Cl_LDL_C")) +
  geom_point(data = combined_data, aes(x = Data_Frame, y = LDL_C, color = "LDL_C")) +
  geom_point(data = combined_data, aes(x = Data_Frame, y = HDL_C, color = "HDL_C")) +
  labs(x = "", y = "Re_C") +
  scale_color_manual(
    values = c("VLDL_C" = "red", "Cl_LDL_C" = "orange", "LDL_C" = "blue", "HDL_C" = "violet") ,
    guide = guide_legend(title = "Cholesterol")) + ggtitle("Correlation of Cholesterol") + theme_bw()

selected_cols <- c("Data_Frame",'Total_C','TL_C_HDL_C','Re_C','VLDL_C','Cl_LDL_C','LDL_C','HDL_C','Total_TG','TG_VLDL', 'TG_LDL','TG_HDL', 'VS_VLDL','C_VS_VLDL','CE_VS_VLDL',
                   'FC_VS_VLDL','TG_VS_VLDL','IDL','CE_VL_HDL','FC_VL_HDL','TG_VL_HDL','L_HDL','C_VL_VLDL','CE_VL_VLDL','FC_VL_VLDL','TG_VL_VLDL','L_VLDL',
                   'C_L_HDL','CE_L_HDL','FC_L_HDL','TG_L_HDL','M_HDL','C_M_LDL','CE_M_LDL','FC_M_LDL','TG_M_LDL','S_LDL','C_S_LDL','CE_S_LDL','FC_S_LDL','TG_S_LDL','VL_HDL','C_VL_HDL')

# Filter only the rows with the selected columns

selected_data <- combined_data %>% select(all_of(selected_cols))

# Calculate the correlation matrix

correlation_matrix <- cor(selected_data[, -1])

# Create a data frame from the correlation matrix

correlation_df <- as.data.frame(as.table(correlation_matrix))

# Rename the columns

colnames(correlation_df) <- c("Variable1", "Variable2", "R2")

# Remove rows where variables are the same (diagonal elements) and remove Data_Frame column

correlation_df <- correlation_df %>%
  filter(Variable1 != Variable2)

# Create a data table of R2 values comparing 'VLDL_C', 'Cl_LDL_C', 'LDL_C', 'HDL_C', and 'C_in_IDL' against each other

R2_data_table <- correlation_df

# Print the correlation matrix and the data table

print("Correlation Matrix:")
print(correlation_matrix)
print("Data Table of R2 Values:")
print(R2_data_table)
correlation_matrix_dt <- datatable(correlation_matrix)

# Create the datatable for the R2 data table

R2_data_table_dt <- datatable(R2_data_table)

# Print the data tables using DT::renderDataTable

correlation_matrix_dt
R2_data_table_dt

# Create the datatable for the R2 data table

R2_data_table_dt <- datatable(R2_data_table)

# Save the R2 data table to a CSV file

write.csv(R2_data_table_dt$x$data, file = "C:/Users/User/Desktop/PhD Documentation/My drafts/Matrix_correlation.csv", row.names = FALSE)