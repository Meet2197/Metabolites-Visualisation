library(randomForest)
library(survival)
library(forestplot)
library(survminer)
library(grid)

# Create labels for forest plot
labels <- c("liverdisease","alcoholicliver","toxicliver","Liverfailure","ChronHepatitis","fibrosecirrho","inflammliver","otherliverD","UnclassLiverd")

# NAs to replace :

all_hazard$liverdisease <-replace_na(all_hazard$liverdisease, 0)
all_hazard$toxicliver <-replace_na(all_hazard$toxicliver, 0)
all_hazard$ChronHepatitis <-replace_na(all_hazard$ChronHepatitis, 0)
all_hazard$Liverfailure <-replace_na(all_hazard$Liverfailure, 0)
all_hazard$fibrosecirrho <-replace_na(all_hazard$fibrosecirrho, 0)
all_hazard$inflammliver <-replace_na(all_hazard$inflammliver, 0)
all_hazard$otherliverD <-replace_na(all_hazard$otherliverD, 0)
all_hazard$UnclassLiverd <-replace_na(all_hazard$UnclassLiverd, 0)
all_hazard$MASLD <-replace_na(all_hazard$MASLD, 0)
all_hazard$MASH <-replace_na(all_hazard$MASH, 0)
all_hazard$metformin <-replace_na(all_hazard$metformin, 0)
all_hazard$alcoholicliver <-replace_na(all_hazard$alcoholicliver, 0)
all_hazard$Diabetes <-replace_na(all_hazard$Diabetes, 0)

# Cox proportional hazard ratio for metformin:

hesin2 <- hesin %>% select('eid','epistart')
hesin3 <- merge(hesin2, init_exam, by = "eid")
hesin3$epidiff <- as.numeric(difftime(hesin3$init_exam, hesin3$epistart, units = "days")) / 365.25
hesin3$epidiff <- abs(hesin3$epidiff)
all_hazard <- merge(hesin3, ALL,  by.x="eid", all.x = TRUE)
all_hazard <- as.data.frame(all_hazard[!duplicated(all_hazard$eid), ])
all_hazard <- merge(all_hazard, liverdisease, by.x="eid", all.x = TRUE)

# Fit Cox proportional hazards model

cox_model1 <- coxph(Surv(epidiff, liverdisease) ~ metformin, data = all_hazard)
cox_model2 <- coxph(Surv(epidiff, alcoholicliver) ~ metformin, data = all_hazard)
cox_model3 <- coxph(Surv(epidiff, toxicliver) ~ metformin, data = all_hazard)
cox_model4 <- coxph(Surv(epidiff, Liverfailure) ~ metformin, data = all_hazard)
cox_model5 <- coxph(Surv(epidiff, ChronHepatitis) ~ metformin, data = all_hazard)
cox_model6 <- coxph(Surv(epidiff, fibrosecirrho) ~ metformin, data = all_hazard)
cox_model7 <- coxph(Surv(epidiff, inflammliver) ~ metformin, data = all_hazard)
cox_model8 <- coxph(Surv(epidiff, otherliverD) ~ metformin, data = all_hazard)
cox_model9 <- coxph(Surv(epidiff, UnclassLiverd) ~ metformin, data = all_hazard)
cox_model10 <- coxph(Surv(epidiff, MASLD) ~ metformin, data = all_hazard)
cox_model11 <- coxph(Surv(epidiff, MRI) ~ metformin, data = all_hazard)

# Extract hazard ratio

hazard_ratio1 <- c(exp(coef(cox_model1)))
hazard_ratio2 <- c(exp(coef(cox_model2)))
hazard_ratio3 <- c(exp(coef(cox_model3)))
hazard_ratio4 <- c(exp(coef(cox_model4)))
hazard_ratio5 <- c(exp(coef(cox_model5)))
hazard_ratio6 <- c(exp(coef(cox_model6)))
hazard_ratio7 <- c(exp(coef(cox_model7)))
hazard_ratio8 <- c(exp(coef(cox_model8)))
hazard_ratio9 <- c(exp(coef(cox_model9)))
hazard_ratio10 <- c(exp(coef(cox_model10)))
hazard_ratio11 <- c(exp(coef(cox_model11)))

# Calculate confidence intervals

conf_interval1 <- rbind(exp(confint(cox_model1)))
conf_interval2 <- rbind(exp(confint(cox_model2)))
conf_interval3 <- rbind(exp(confint(cox_model3)))
conf_interval4 <- rbind(exp(confint(cox_model4)))
conf_interval5 <- rbind(exp(confint(cox_model5)))
conf_interval6 <- rbind(exp(confint(cox_model6))) 
conf_interval7 <- rbind(exp(confint(cox_model7))) 
conf_interval8 <- rbind(exp(confint(cox_model8)))
conf_interval9 <- rbind(exp(confint(cox_model9)))
conf_interval10 <- rbind(exp(confint(cox_model10)))
conf_interval11 <- rbind(exp(confint(cox_model11)))

# Extract p-values

p_value1 <- c(summary(cox_model1)$coefficients["metformin", "Pr(>|z|)"])
p_value2 <- c(summary(cox_model2)$coefficients["metformin", "Pr(>|z|)"])
p_value3 <- c(summary(cox_model3)$coefficients["metformin", "Pr(>|z|)"])
p_value4 <- c(summary(cox_model4)$coefficients["metformin", "Pr(>|z|)"])
p_value5 <- c(summary(cox_model5)$coefficients["metformin", "Pr(>|z|)"])
p_value6 <- c(summary(cox_model6)$coefficients["metformin", "Pr(>|z|)"])
p_value7 <- c(summary(cox_model7)$coefficients["metformin", "Pr(>|z|)"])
p_value8 <- c(summary(cox_model8)$coefficients["metformin", "Pr(>|z|)"])
p_value9 <- c(summary(cox_model9)$coefficients["metformin", "Pr(>|z|)"])
p_value10 <- c(summary(cox_model10)$coefficients["metformin", "Pr(>|z|)"])
p_value11 <- c(summary(cox_model11)$coefficients["metformin", "Pr(>|z|)"])

# Cox proportional hazard ratio for Pioglitazone:

all_hazard2 <- merge(hesin3, ALL_pioglitazone,  by.x="eid", all.x = TRUE)
all_hazard2 <- as.data.frame(all_hazard2[!duplicated(all_hazard2$eid), ])

# Fit Cox proportional hazards model

cox_model_pio1 <- coxph(Surv(epidiff, liverdisease) ~ pioglitazone, data = all_hazard2)
cox_model_pio2 <- coxph(Surv(epidiff, alcoholicliver) ~ pioglitazone, data = all_hazard2)
cox_model_pio3 <- coxph(Surv(epidiff, toxicliver) ~ pioglitazone, data = all_hazard2)
cox_model_pio4 <- coxph(Surv(epidiff, Liverfailure) ~ pioglitazone, data = all_hazard2)
cox_model_pio5 <- coxph(Surv(epidiff, ChronHepatitis) ~ pioglitazone, data = all_hazard2)
cox_model_pio6 <- coxph(Surv(epidiff, fibrosecirrho) ~ pioglitazone, data = all_hazard2)
cox_model_pio7 <- coxph(Surv(epidiff, inflammliver) ~ pioglitazone, data = all_hazard2)
cox_model_pio8 <- coxph(Surv(epidiff, otherliverD) ~ pioglitazone, data = all_hazard2)
cox_model_pio9 <- coxph(Surv(epidiff, UnclassLiverd) ~ pioglitazone, data = all_hazard2)
cox_model_pio10 <- coxph(Surv(epidiff, MASLD) ~ pioglitazone, data = all_hazard2)
cox_model_pio11 <- coxph(Surv(epidiff, MRI) ~ pioglitazone, data = all_hazard2)


# Extract hazard ratio

hazards_ratio1 <- c(exp(coef(cox_model_pio1)))
hazards_ratio2 <- c(exp(coef(cox_model_pio2)))
hazards_ratio3 <- c(exp(coef(cox_model_pio3)))
hazards_ratio4 <- c(exp(coef(cox_model_pio4)))
hazards_ratio5 <- c(exp(coef(cox_model_pio5)))
hazards_ratio6 <- c(exp(coef(cox_model_pio6)))
hazards_ratio7 <- c(exp(coef(cox_model_pio7))) 
hazards_ratio8 <- c(exp(coef(cox_model_pio8)))
hazards_ratio9 <- c(exp(coef(cox_model_pio9)))
hazards_ratio10 <- c(exp(coef(cox_model_pio10)))
hazards_ratio11 <- c(exp(coef(cox_model_pio11)))


# Calculate confidence intervals
conf_intervals1 <- rbind(exp(confint(cox_model_pio1)))
conf_intervals2 <- rbind(exp(confint(cox_model_pio2)))
conf_intervals3 <- rbind(exp(confint(cox_model_pio3)))
conf_intervals4 <- rbind(exp(confint(cox_model_pio4)))
conf_intervals5 <- rbind(exp(confint(cox_model_pio5)))
conf_intervals6 <- rbind(exp(confint(cox_model_pio6)))
conf_intervals7 <- rbind(exp(confint(cox_model_pio7)))
conf_intervals8 <- rbind(exp(confint(cox_model_pio8)))
conf_intervals9 <- rbind(exp(confint(cox_model_pio9)))
conf_intervals10 <- rbind(exp(confint(cox_model_pio10)))
conf_intervals11 <- rbind(exp(confint(cox_model_pio11)))

# Extract p-values
p_values1 <- c(summary(cox_model_pio1)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values2 <- c(summary(cox_model_pio2)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values3 <- c(summary(cox_model_pio3)$coefficients["pioglitazone", "Pr(>|z|)"]) 
p_values4 <- c(summary(cox_model_pio4)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values5 <- c(summary(cox_model_pio5)$coefficients["pioglitazone", "Pr(>|z|)"]) 
p_values6 <- c(summary(cox_model_pio6)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values7 <- c(summary(cox_model_pio7)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values8 <- c(summary(cox_model_pio8)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values9 <- c(summary(cox_model_pio9)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values10 <- c(summary(cox_model_pio10)$coefficients["pioglitazone", "Pr(>|z|)"])
p_values11 <- c(summary(cox_model_pio11)$coefficients["pioglitazone", "Pr(>|z|)"])

# Cox proportional hazard ratio for Pioglitazone:

all_hazard3 <- merge(hesin3, ALL_ramipril,  by.x="eid", all.x = TRUE)
all_hazard3 <- as.data.frame(all_hazard3[!duplicated(all_hazard3$eid), ])

# Fit Cox proportional hazards model

cox_model_ram1 <- coxph(Surv(epidiff, liverdisease) ~ ramipril, data = all_hazard3)
cox_model_ram2 <- coxph(Surv(epidiff, alcoholicliver) ~ ramipril, data = all_hazard3)
cox_model_ram3 <- coxph(Surv(epidiff, toxicliver) ~ ramipril, data = all_hazard3)
cox_model_ram4 <- coxph(Surv(epidiff, Liverfailure) ~ ramipril, data = all_hazard3)
cox_model_ram5 <- coxph(Surv(epidiff, ChronHepatitis) ~ ramipril, data = all_hazard3)
cox_model_ram6 <- coxph(Surv(epidiff, fibrosecirrho) ~ ramipril, data = all_hazard3)
cox_model_ram7 <- coxph(Surv(epidiff, inflammliver) ~ ramipril, data = all_hazard3)
cox_model_ram8 <- coxph(Surv(epidiff, otherliverD) ~ ramipril, data = all_hazard3)
cox_model_ram9 <- coxph(Surv(epidiff, UnclassLiverd) ~ ramipril, data = all_hazard3)

# Extract hazard ratio

hazards_ratio1 <- c(exp(coef(cox_model_ram1)))
hazards_ratio2 <- c(exp(coef(cox_model_ram2)))
hazards_ratio3 <- c(exp(coef(cox_model_ram3)))
hazards_ratio4 <- c(exp(coef(cox_model_ram4)))
hazards_ratio5 <- c(exp(coef(cox_model_ram5)))
hazards_ratio6 <- c(exp(coef(cox_model_ram6)))
hazards_ratio7 <- c(exp(coef(cox_model_ram7))) 
hazards_ratio8 <- c(exp(coef(cox_model_ram8)))
hazards_ratio9 <- c(exp(coef(cox_model_ram9)))

# Calculate confidence intervals
confi_intervals1 <- rbind(exp(confint(cox_model_ram1)))
confi_intervals2 <- rbind(exp(confint(cox_model_ram2)))
confi_intervals3 <- rbind(exp(confint(cox_model_ram3)))
confi_intervals4 <- rbind(exp(confint(cox_model_ram4)))
confi_intervals5 <- rbind(exp(confint(cox_model_ram5)))
confi_intervals6 <- rbind(exp(confint(cox_model_ram6)))
confi_intervals7 <- rbind(exp(confint(cox_model_ram7)))
confi_intervals8 <- rbind(exp(confint(cox_model_ram8)))
confi_intervals9 <- rbind(exp(confint(cox_model_ram9)))

# Extract p-values
P_valuesS1 <- c(summary(cox_model_ram1)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS2 <- c(summary(cox_model_ram2)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS3 <- c(summary(cox_model_ram3)$coefficients["ramipril", "Pr(>|z|)"]) 
P_valuesS4 <- c(summary(cox_model_ram4)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS5 <- c(summary(cox_model_ram5)$coefficients["ramipril", "Pr(>|z|)"]) 
P_valuesS6 <- c(summary(cox_model_ram6)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS7 <- c(summary(cox_model_ram7)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS8 <- c(summary(cox_model_ram8)$coefficients["ramipril", "Pr(>|z|)"])
P_valuesS9 <- c(summary(cox_model_ram9)$coefficients["ramipril", "Pr(>|z|)"])

# Repeat for other models...

surv_fit1 <- survfit(cox_model2)
surv_fit2 <- survfit(cox_model3)
summary_surv_fit1 <- summary(surv_fit1)# Create summary of survival fits
summary_surv_fit2 <- summary(surv_fit2)
time_shift <- min(summary_surv_fit1$time, summary_surv_fit2$time) # Assuming surv_fit1 and surv_fit2 are already defined
surv_prob1 <- summary_surv_fit1$surv # Create survival curves
time_points1 <- summary_surv_fit1$time - time_shift
surv_prob2 <- summary_surv_fit2$surv
time_points2 <- summary_surv_fit2$time - time_shift

# Plot survival curves

ggplot() + 
  geom_step(data = data.frame(time = time_points1, surv = surv_prob1, group = "alcoholicliver"), aes(x = time, y = surv, group = group), color = "red", linewidth = 1) +
  geom_step(data = data.frame(time = time_points2, surv = surv_prob2, group = "Toxicliver"), aes(x = time, y = surv, group = group), color = "blue", linewidth = 1) +
  labs(x = "Time(Years)", y = "Survival Probability", title = "Survival Curves for Cox Models", color = "Group") +  
  theme_minimal()   # Add more geom_step() for other survival curves...

# forest plot
average_mean <- mean(hazard_ratios)
forest_data <- forest_data %>%
  mutate(mean_limited = ifelse(mean > 3, 3, mean),
         upper_limited = ifelse(mean > 3, 3, upper),
         arrow_label = ifelse(mean > 3, " > 3", ""))  # Add arrow label for hazard ratios > 3

# Plot forest plot using ggplot2 with modifications
arrow_data <- forest_data %>%
  filter(upper > 3)

# Plot forest plot using ggplot2 with modifications
ggplot(forest_data, aes(x = mean_limited, y = label)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), color = "red") +  # Add red color to confidence interval lines
  geom_vline(xintercept = average_mean / 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Line of null effect
  geom_text(aes(label = arrow_label), hjust = -0.5, color = "red", size = 5, check_overlap = TRUE, fontface = "bold") +  # Add arrow label with red bold font
  geom_segment(data = arrow_data, aes(x = upper_limited + 0.1, xend = upper_limited + 0.5, y = label, yend = label),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "black", size = 0.5) +  # Add arrow segments
  labs(x = "Hazard Ratio", title = "Forest Plot of Hazard Ratio for metformin consumption") +  # Add title
  theme_minimal() +
  theme(panel.grid = element_blank())   # Remove grid lines

# Combine hazard ratios, confidence intervals, and p-values into a data frame
forest_data2 <- data.frame(
  treatment = c("liverdisease","alcoholicliver","toxicliver","Liverfailure","ChronHepatitis","fibrosecirrho","inflammliver","otherliverD","UnclassLiverd","MASLD","MRI"),
  mean = c(hazards_ratio1, hazards_ratio2, hazards_ratio3, hazards_ratio4, hazards_ratio5, hazards_ratio6, hazards_ratio7, hazards_ratio8, hazards_ratio9, hazards_ratio10, hazards_ratio11),
  lower = c(conf_intervals1[1], conf_intervals2[1], conf_intervals3[1], conf_intervals4[1], conf_intervals5[1], conf_intervals6[1], conf_intervals7[1], conf_intervals8[1], conf_intervals9[1], conf_intervals10[1], conf_intervals11[1]),
  upper = c(conf_intervals1[2], conf_intervals2[2], conf_intervals3[2], conf_intervals4[2], conf_intervals5[2], conf_intervals6[2], conf_intervals7[2], conf_intervals8[2], conf_intervals9[2], conf_intervals10[2],conf_intervals11[2]),
  p_value = c(p_values1, p_values2, p_values3, p_values4, p_values5, p_values6, p_values7, p_values8, p_values9, p_values10, p_values11))

# Calculate average mean for reference
average_mean2 <- mean(forest_data2$mean)

# Plot forest plot using ggplot2 with modifications
ggplot(forest_data2, aes(x = mean, y = treatment)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), color = "red") +  # Add red color to confidence interval lines
  geom_vline(xintercept = average_mean2 / 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Line of null effect
  geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), hjust = -0.5, color = "black", size = 3.5, check_overlap = TRUE) +  # Add asterisk for significant p-values
  labs(x = "Hazard Ratio", y = "Treatment", title = "Forest Plot of Pioglitazone Hazard Ratio") +  # Add title and axis labels
  theme_minimal() +
  theme(panel.grid = element_blank())
# 
# Kaplan Meier: 
# Fit the Cox proportional hazards model:
cox_model_sur <- coxph(Surv(epidiff, liverdisease) ~ death, data = all_hazard)

# Create Kaplan-Meier curve
surv_object <- Surv(time = all_hazard$epidiff, event = all_hazard$liverdisease)
km_fit <- survfit(surv_object ~ death, data = all_hazard)

# Plot the Kaplan-Meier curve with a legend
ggsurvplot(km_fit, data = all_hazard, risk.table = TRUE, ggtheme = theme_minimal(), legend.labs = c("Survived", "death"))
