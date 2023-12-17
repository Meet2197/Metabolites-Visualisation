# 1.1 Extract specific variables from each dataframe (LDL Comparison)

S_VLDL <- healthy_df$S_VLDL
L_LDL <- metformin_filtered_df$L_LDL
M_LDL <- NASH_1$M_LDL

# Create a combined dataframe of Metabolites with "Healthy", "Metformin Filtered", "NASH", "MAFLD" datasets.  
# combined data = "Healthy", "Metformin Filtered", "NASH", "MAFLD" dataframes mentioned

comparison_data20 <- data.frame(
  Variables = rep(c("S_VLDL", "L_LDL", "M_LDL"), each = 4),
  Category = rep(c("Healthy", "Metformin Filtered", "NASH", "MAFLD"), times = 3),
  Mean = c(mean(S_VLDL), mean(L_LDL), mean(M_LDL)),
  SD = c(sd(S_VLDL), sd(L_LDL), sd(M_LDL)),
  Median = c(median(S_VLDL), median(L_LDL), median(M_LDL)),
  Min = c(min(S_VLDL), min(L_LDL), min(M_LDL)),
  Max = c(max(S_VLDL), max(L_LDL), max(M_LDL)),
  IQR = c(IQR(S_VLDL), IQR(L_LDL), IQR(M_LDL)),
  Q1 = c(quantile(S_VLDL, 0.25), quantile(L_LDL, 0.25), quantile(M_LDL, 0.25)),
  Q3 = c(quantile(S_VLDL, 0.75), quantile(L_LDL, 0.75), quantile(M_LDL, 0.75))
)

# Plotting the whisker box plot of comparing one major metabolite with other sub-metabolites 

ggplot(combined_data, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("LDL Comparison(S_VLDL,L_LDL,M_LDL)") +
  scale_fill_manual(values = c("S_VLDL" = "lightblue", "L_LDL" = "lightgreen", "M_LDL" = "lightyellow"))

# Co-relation of two different Statstical overview of one metabolite

correlation1 <- cor(combined_data$Mean, combined_data$Median)
correlation2 <- cor(combined_data$SD, combined_data$IQR)
correlation3 <- cor(combined_data$Min, combined_data$Max)

# Perform linear regression of Combined dataset of 

regression1 <- lm(Mean ~ Median, data = combined_data)
regression2 <- lm(SD ~ IQR, data = combined_data)
regression3 <- lm(Min ~ Max, data = combined_data)

# ggplot combined data 

ggplot(combined_data, aes(x = Median, y = Mean)) +
  geom_point(aes(fill = Variables), shape = 23, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "orange", formula = y ~ x, data = combined_data) +
  labs(x = "Median", y = "Mean") +
  ggtitle("Linear Regression Scatter Plot") +
  scale_fill_manual(values = c("S_VLDL" = "lightblue", "L_LDL" = "lightgreen", "M_LDL" = "lightyellow")) +
  annotate("text", x = max(combined_data$Median), y = max(combined_data$Mean), label = paste("Correlation 1:", round(correlation1, 2)), hjust = 1, vjust = -0.5, color = "black") +
  annotate("text", x = max(combined_data$Median), y = max(combined_data$Mean), label = paste("Correlation 2:", round(correlation2, 2)), hjust = 1, vjust = -1, color = "black") +
  annotate("text", x = max(combined_data$Median), y = max(combined_data$Mean), label = paste("Correlation 3:", round(correlation3, 2)), hjust = 1, vjust = -1.5, color = "black")


# Plotting the whisker box plot with correlation and regression lines

ggplot(combined_data, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "orange", formula = y ~ x, data = combined_data) +
  labs(x = "Category", y = "Mean") +
  ggtitle("LDL Comparison (S_VLDL, L_LDL, M_LDL)") +
  scale_fill_manual(values = c("S_VLDL" = "lightblue", "L_LDL" = "lightgreen", "M_LDL" = "lightyellow")) +
  annotate("text", x = "NASH", y = max(combined_data$Mean), label = paste("Correlation1:", round(correlation1, 2)), hjust = 0, vjust = -0.5, color = "black") +
  annotate("text", x = "NASH", y = max(combined_data$Mean) - 1, label = paste("Correlation2:", round(correlation2, 2)), hjust = 0, vjust = -0.5, color = "black") +
  annotate("text", x = "NASH", y = max(combined_data$Mean) - 2, label = paste("Correlation3:", round(correlation3, 2)), hjust = 0, vjust = -0.5, color = "black") +
  geom_abline(intercept = coef(regression1)[1], slope = coef(regression1)[2], color = "red", linetype = "dashed") +
  geom_abline(intercept = coef(regression2)[1], slope = coef(regression2)[2], color = "blue", linetype = "dashed") +
  geom_abline(intercept = coef(regression3)[1], slope = coef(regression3)[2], color = "green", linetype = "dashed")


# 2 Extract specific variables from each dataframe (LDL Comparison)

LDL_C <- healthy_df$LDL_C
VLDL_C <- metformin_filtered_df$VLDL_C
Re_C <- NASH_1$Re_C
Tl_C <- MAFLD$Tl_C

# Create a combined dataframe with variables to compare 

comparison_data3 <- data.frame(
  Variables = rep(c("LDL_C", "VLDL_C", "Re_C", "Tl_C"), each = 4),
  Category = rep(c("Healthy", "Metformin Filtered", "NASH", "MAFLD"), times = 4),
  Mean = c(mean(LDL_C), mean(VLDL_C), mean(Re_C), mean(Tl_C)),
  SD = c(sd(LDL_C), sd(VLDL_C), sd(Re_C), sd(Tl_C)),
  Median = c(median(LDL_C), median(VLDL_C), median(Re_C), median(Tl_C)),
  Min = c(min(LDL_C), min(VLDL_C), min(Re_C), min(Tl_C)),
  Max = c(max(LDL_C), max(VLDL_C), max(Re_C), max(Tl_C)),
  IQR = c(IQR(LDL_C), IQR(VLDL_C), IQR(Re_C), IQR(Tl_C)),
  Q1 = c(quantile(LDL_C, 0.25), quantile(VLDL_C, 0.25), quantile(Re_C, 0.25), quantile(Tl_C, 0.25)),
  Q3 = c(quantile(LDL_C, 0.75), quantile(VLDL_C, 0.75), quantile(Re_C, 0.75), quantile(Tl_C, 0.75))
)

# Plotting the whisker box plot of dataframes

ggplot(comparison_data3, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("Cholesterol Comparison (LDL_C, VLDL_C, Re_C, Tl_C)") +
  scale_fill_manual(values = c("LDL_C" = "lightblue", "VLDL_C" = "lightgreen", "Re_C" = "lightyellow", "Tl_C" = "lightpink"))

correlation <- cor(comparison_data3$Mean, comparison_data3$Median)

# Perform linear regression
regression <- lm(Mean ~ Median, data = comparison_data3)

# Plotting the whisker box plot with correlation and regression lines
ggplot(comparison_data3, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "orange", formula = y ~ x, data = comparison_data3) +
  labs(x = "Category", y = "Mean") +
  ggtitle("Cholesterol Comparison (LDL_C, VLDL_C, Re_C, Tl_C)") +
  scale_fill_manual(values = c("LDL_C" = "lightblue", "VLDL_C" = "lightgreen", "Re_C" = "lightyellow", "Tl_C" = "lightpink")) +
  annotate("text", x = "NASH", y = max(comparison_data3$Mean), label = paste("Correlation:", round(correlation, 2)), hjust = 0, vjust = -0.5, color = "black")

# 3 Extract specific variables from each dataframe (LDL Comparison)

Apo_B <- healthy_df$Apo_B
Apo_B <- metformin_filtered_df$Apo_B

# Create a combined dataframe
comparison_data <- data.frame(
  Variables = rep(c("Apo_B", "Apo_A1"), each = 4),
  Category = rep(c("Healthy", "Metformin Filtered", "NASH", "MAFLD"), times = 2),
  Mean = c(mean(Apo_B), mean(Apo_A1)),
  SD = c(sd(Apo_B), sd(Apo_A1)),
  Median = c(median(Apo_B), median(Apo_A1)),
  Min = c(min(Apo_B), min(Apo_A1)),
  Max = c(max(Apo_B), max(Apo_A1)),
  IQR = c(IQR(Apo_B), IQR(Apo_A1)),
  Q1 = c(quantile(Apo_B, 0.25), quantile(Apo_A1, 0.25)),
  Q3 = c(quantile(Apo_B, 0.75), quantile(Apo_A1, 0.75))
)

ggplot(comparison_data, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("Enzyme Comparison (Apo_B, Apo_A1)") +
  scale_fill_manual(values = c("Apo_B" = "lightblue", "Apo_A1" = "lightgreen"))


correlation <- cor(comparison_data1$Mean, comparison_data1$Median)

# Perform linear regression
regression <- lm(Mean ~ Median, data = comparison_data1)

# Plotting the whisker box plot with correlation and regression lines

ggplot(comparison_data1, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "orange", formula = y ~ x, data = comparison_data1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("Enyzme Comparison (Apo_B, Apo_A1)") +
  scale_fill_manual(values = c("Apo_B" = "lightblue", "Apo_A1" = "lightgreen")) +
  annotate("text", x = "NASH", y = max(comparison_data1$Mean), label = paste("Correlation:", round(correlation, 2)), hjust = 0, vjust = -0.5, color = "black")

# 4 Extract specific variables from each dataframe (LDL Comparison)

IDL <- healthy_df$IDL
C_in_IDL <- metformin_filtered_df$C_in_IDL
CE_in_IDL <- MAFLD$CE_in_IDL

# Create a combined dataframe
comparison_data6 <- data.frame(
  Variables = rep(c("IDL", "C_in_IDL", "CE_in_IDL"), each = 4),
  Category = rep(c("Healthy", "Metformin Filtered", "NASH", "MAFLD"), times = 3),
  Mean = c(mean(IDL), mean(C_in_IDL), mean(CE_in_IDL)),
  SD = c(sd(IDL), sd(C_in_IDL), sd(CE_in_IDL)),
  Median = c(median(IDL), median(C_in_IDL), median(CE_in_IDL)),
  Min = c(min(IDL), min(C_in_IDL), min(CE_in_IDL)),
  Max = c(max(IDL), max(C_in_IDL), max(CE_in_IDL)),
  IQR = c(IQR(IDL), IQR(C_in_IDL), IQR(CE_in_IDL)),
  Q1 = c(quantile(IDL, 0.25), quantile(C_in_IDL, 0.25), quantile(CE_in_IDL, 0.25)),
  Q3 = c(quantile(IDL, 0.75), quantile(C_in_IDL, 0.75), quantile(CE_in_IDL, 0.75))
)

ggplot(comparison_data6, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("IDL Comparison (IDL, C_in_IDL, CE_in_IDL)") +
  scale_fill_manual(values = c("IDL" = "lightblue", "C_in_IDL" = "lightgreen", "CE_in_IDL" = "lightyellow"))


# 5 Extract specific variables from each dataframe (LDL Comparison)

DHA <- healthy_df$DHA
Glutamine <- metformin_filtered_df$Glutamine
Glycine <- MAFLD$Glycine
Glucose <- Nash$Glucose

# Create a combined dataframe
comparison_data7 <- data.frame(
  Variables = rep(c("DHA", "Glutamine", "Glycine", "Glucose"), each = 4),
  Category = rep(c("Healthy", "Metformin Filtered", "MAFLD", "NASH"), times = 4),
  Mean = c(mean(DHA), mean(Glutamine), mean(Glycine), mean(Glucose)),
  SD = c(sd(DHA), sd(Glutamine), sd(Glycine), sd(Glucose)),
  Median = c(median(DHA), median(Glutamine), median(Glycine), median(Glucose)),
  Min = c(min(DHA), min(Glutamine), min(Glycine), min(Glucose)),
  Max = c(max(DHA), max(Glutamine), max(Glycine), max(Glucose)),
  IQR = c(IQR(DHA), IQR(Glutamine), IQR(Glycine), IQR(Glucose)),
  Q1 = c(quantile(DHA, 0.25), quantile(Glutamine, 0.25), quantile(Glycine, 0.25), quantile(Glucose, 0.25)),
  Q3 = c(quantile(DHA, 0.75), quantile(Glutamine, 0.75), quantile(Glycine, 0.75), quantile(Glucose, 0.75))
)

# Plotting the whisker box plot
ggplot(comparison_data7, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("AAs Comparison (DHA, Glutamine, Glycine, Glucose)") +
  scale_fill_manual(values = c("DHA" = "lightblue", "Glutamine" = "lightgreen", "Glycine" = "lightyellow", "Glucose" = "lightpink"))

correlation1 <- cor(combined_data$Mean, combined_data$Median)
correlation2 <- cor(combined_data$SD, combined_data$IQR)
correlation3 <- cor(combined_data$Min, combined_data$Max)

# Perform linear regression
regression1 <- lm(Mean ~ Median, data = combined_data)
regression2 <- lm(SD ~ IQR, data = combined_data)
regression3 <- lm(Min ~ Max, data = combined_data)

# Plotting the whisker box plot with correlation and regression lines
ggplot(comparison_data1, aes(x = Category, y = Mean, fill = Variables)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, color = "red") +
  geom_point(aes(y = Median), shape = 23, size = 3, fill = "black") +
  geom_segment(aes(y = Q1, yend = Q3, xend = Category), color = "blue", size = 1) +
  geom_segment(aes(y = Mean - SD, yend = Mean + SD, xend = Category), color = "purple", size = 1) +
  geom_segment(aes(y = Mean - IQR/2, yend = Mean + IQR/2, xend = Category), color = "green", size = 1) +
  labs(x = "Category", y = "Mean") +
  ggtitle("CE Comparison (DHA, Glutamine, Glycine, Glucose)") +
  scale_fill_manual(values = c("DHA" = "lightblue", "Glutamine" = "lightgreen", "Glycine" = "lightyellow", "Glucose" = "lightpink")) +
  annotate("text", x = "NASH", y = max(comparison_data1$Mean), label = paste("Correlation:", round(correlation, 2)), hjust = 0, vjust = -0.5, color = "black")
