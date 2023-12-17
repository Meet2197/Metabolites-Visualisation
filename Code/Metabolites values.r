
# Calculate summary statistics for each metabolites e.g VLDL_C , 

Cl_LDL_C_summary <- data.frame(
  Variables = "Cl_LDL_C",
  Mean = round(mean(Cl_LDL_C, na.rm = TRUE), 2),
  SD = round(sd(Cl_LDL_C, na.rm = TRUE), 2),
  Median = round(median(Cl_LDL_C, na.rm = TRUE), 2),
  Q1 = round(quantile(Cl_LDL_C, probs = 0.25, na.rm = TRUE), 2),
  Q3 = round(quantile(Cl_LDL_C, probs = 0.75, na.rm = TRUE), 2),
  IQR = round(IQR(Cl_LDL_C, na.rm = TRUE), 2),
  Min = round(min(Cl_LDL_C, na.rm = TRUE), 2),
  Max = round(max(Cl_LDL_C, na.rm = TRUE), 2))

dataframes <- list(metformin_df, healthy_df, no_liver_disease, Fibrosis_cirrhosis, MAFLD_filtered)

for (df in dataframes) {
  summary_row <- data.frame(
    Variables = "Tl_C",
    Mean = round(mean(df$Tl_C, na.rm = TRUE), 2),
    SD = round(sd(df$Tl_C, na.rm = TRUE), 2),
    Median = round(median(df$Tl_C, na.rm = TRUE), 2),
    Q1 = round(quantile(df$Tl_C, probs = 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(df$Tl_C, probs = 0.75, na.rm = TRUE), 2),
    IQR = round(IQR(df$Tl_C, na.rm = TRUE), 2),
    Min = round(min(df$Tl_C, na.rm = TRUE), 2),
    Max = round(max(df$Tl_C, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
  
  Tl_C_summary <- rbind(Tl_C_summary, summary_row)
}

print(Tl_C_summary)


Re_C_summary <- data.frame(
  Variables = "Re_C",
  Mean = round(mean(Re_C, na.rm = TRUE), 2),
  SD = round(sd(Re_C, na.rm = TRUE), 2),
  Median = round(median(Re_C, na.rm = TRUE), 2),
  Q1 = round(quantile(Re_C, probs = 0.25, na.rm = TRUE), 2),
  Q3 = round(quantile(Re_C, probs = 0.75, na.rm = TRUE), 2),
  IQR = round(IQR(Re_C, na.rm = TRUE), 2),
  Min = round(min(Re_C, na.rm = TRUE), 2),
  Max = round(max(Re_C, na.rm = TRUE), 2))

# Plot boxplot using boxplot parameters

Cl_LDL_C$Variables <- as.factor(vldl_c_summary$Variables)
Cl_LDL_C_summary %>% 
  pivot_longer(cols = -Variables, names_to = "Statistic", values_to = "Value") %>%
  ggplot(aes(x = Statistic, y = Value, fill = Variables)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Cl_LDL_C") +
  ggtitle("Cl_LDL_C")

