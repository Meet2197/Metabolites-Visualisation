library(randomForest)
library(survival)

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

# Extract hazard ratio

hazard_ratios <- c(exp(coef(cox_model1)),exp(coef(cox_model2)), exp(coef(cox_model3)),exp(coef(cox_model4)), exp(coef(cox_model5)),
                   exp(coef(cox_model6)), exp(coef(cox_model7)), exp(coef(cox_model8)), exp(coef(cox_model9)))

# Calculate confidence intervals
conf_intervals <- rbind(exp(confint(cox_model1)), exp(confint(cox_model2)), exp(confint(cox_model3)), exp(confint(cox_model4)),
                        exp(confint(cox_model5)), exp(confint(cox_model6)), exp(confint(cox_model7)), exp(confint(cox_model8)),exp(confint(cox_model9)))

# Extract p-values

p_values <- c( summary(cox_model1)$coefficients["metformin", "Pr(>|z|)"], summary(cox_model2)$coefficients["metformin", "Pr(>|z|)"],
               summary(cox_model3)$coefficients["metformin", "Pr(>|z|)"], summary(cox_model4)$coefficients["metformin", "Pr(>|z|)"],
               summary(cox_model5)$coefficients["metformin", "Pr(>|z|)"], summary(cox_model6)$coefficients["metformin", "Pr(>|z|)"],
               summary(cox_model7)$coefficients["metformin", "Pr(>|z|)"], summary(cox_model8)$coefficients["metformin", "Pr(>|z|)"],
               summary(cox_model9)$coefficients["metformin", "Pr(>|z|)"])

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

# Extract hazard ratio

hazard_ratios <- c(exp(coef(cox_model_pio1)),exp(coef(cox_model_pio2)), exp(coef(cox_model_pio3)),exp(coef(cox_model_pio4)), exp(coef(cox_model_pio5)),
                   exp(coef(cox_model_pio6)), exp(coef(cox_model_pio7)), exp(coef(cox_model_pio8)), exp(coef(cox_model_pio9)))

# Calculate confidence intervals
conf_intervals <- rbind(exp(confint(cox_model_pio1)), exp(confint(cox_model_pio2)), exp(confint(cox_model_pio3)), exp(confint(cox_model_pio4)),
                        exp(confint(cox_model_pio5)), exp(confint(cox_model_pio6)), exp(confint(cox_model_pio7)), exp(confint(cox_model_pio8)),exp(confint(cox_model_pio9)))

# Extract p-values
p_values <- c( summary(cox_model_pio1)$coefficients["pioglitazone", "Pr(>|z|)"], summary(cox_model_pio2)$coefficients["pioglitazone", "Pr(>|z|)"],
               summary(cox_model_pio3)$coefficients["pioglitazone", "Pr(>|z|)"], summary(cox_model_pio4)$coefficients["pioglitazone", "Pr(>|z|)"],
               summary(cox_model_pio5)$coefficients["pioglitazone", "Pr(>|z|)"], summary(cox_model_pio6)$coefficients["pioglitazone", "Pr(>|z|)"],
               summary(cox_model_pio7)$coefficients["pioglitazone", "Pr(>|z|)"], summary(cox_model_pio8)$coefficients["pioglitazone", "Pr(>|z|)"],
               summary(cox_model_pio9)$coefficients["pioglitazone", "Pr(>|z|)"])


# Fit survival curve
surv_obj <- survfit(Surv(epidiff, death) ~ 1, data = all_hazard2)

# Plot survival curve with customized settings
plot(surv_obj, xlab = "Years", ylab = "Survival Probability", main = "Survival Curve for Patients with Liver Disease",
     col = "blue", # Change line color
     lwd = 2,      # Increase line width
     lty = 1)      # Change line type

# Plot survival curves

ggplot() + 
  geom_step(data = data.frame(time = time_points1, surv = surv_prob1, group = "alcoholicliver"), aes(x = time, y = surv, group = group), color = "red", linewidth = 1) +
  geom_step(data = data.frame(time = time_points2, surv = surv_prob2, group = "Toxicliver"), aes(x = time, y = surv, group = group), color = "blue", linewidth = 1) +
  labs(x = "Time", y = "Survival Probability", title = "Survival Curves for Cox Models", color = "Group") +  
  theme_minimal()   # Add more geom_step() for other survival curves...

# Create forest plot

forestplot( labeltext = labels, mean = hazard_ratios, lower = conf_intervals[, 1], upper = conf_intervals[, 2],
            pvalues = p_values, txt_gp = fpTxtGp(label = gpar(fontsize = 12, col = "blue")),
            col = fpColors(box = "black", lines = "red", text = "black"), xlab = "Hazard Ratio", graphwidth = unit(10, "cm"),
            boxsize = 0.3, meanlines = TRUE, xticks = c(0.1, 0.5, 1, 2, 3), colgap = unit(1, "cm"), lwd.ci = 2, lwd.mean = 3, lwd.p = 2,
            lwd.zero = 1, clip = c(0.1, 3), pval = 0.05, symbol = "circle", new_page = TRUE, vertices = TRUE,vertex.placement = c(0, 1)) 

grid.lines(x = unit(seq(0, 1, length.out = length(labels)), "npc"), y = unit(0, "npc"), gp = gpar(col = "blue", lwd = 1, lty = "dotted"))


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
