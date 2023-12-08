library(patchwork)

# grouping variable

Dataframes <- c(rep("MAFLD", nrow(MAFLD)),
                     rep("Nash", nrow()),
                     rep("Others Medications", nrow(nonmetformin_df)),
                     rep("Metformin (Filtered)", nrow(metformin_df)),
                     rep("Healthy", nrow(Healthy_Patient)))

# 2.1 Cl_LDL_C values from all dataframes
Cl_LDL_C_values <- c(MAFLD$Cl_LDL_C,Nash$Cl_LDL_C, nonmetformin_df$Cl_LDL_C, metformin_filtered_df$Cl_LDL_C, healthy_df$Cl_LDL_C)

# 2.2 dataframe with combined values and grouping variable
combined_data2 <- data.frame(Cl_LDL_C = Cl_LDL_C_values, Dataframe = dataframe_names)

# 2.3 Calculate whisker values

whisker_data1 <- combined_data2 %>%
  group_by(Dataframe) %>%
  summarise(lower_whisker = quantile(Cl_LDL_C, 0.25) - 1.5 * IQR(Cl_LDL_C),
            upper_whisker = quantile(Cl_LDL_C, 0.75) + 1.5 * IQR(Cl_LDL_C))
# 
# 2.4 'whisker_data' to a tibble

whisker_data2 <- as_tibble(whisker_data1)

# 2.5 box plot with whiskers
ggplot() +
  geom_boxplot(data = combined_data2, aes(x = Dataframe, y = Cl_LDL_C)) +
  geom_errorbar(data = whisker_data2, aes(x = Dataframe, ymin = lower_whisker, ymax = upper_whisker), width = 0.2, color = "red") +
  labs(x = "Dataframe", y = "Cl_LDL_C") +
  ggtitle("Clinical_LDL_Cholesterol Values Comparison")
