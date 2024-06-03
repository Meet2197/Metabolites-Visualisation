library(car)

# Model for regression : 

model <- glm(pioglitazone ~ Diabetes + MASLD + MASH + alcoholicliver + toxicliver + Liverfailure + ChronHepatitis + fibrosecirrho + inflammliver + otherliverD + UnclassLiverd, data = pioglitazone_psm1, family = binomial)

# Summarize the model

summary(model)

# Confidence Intervals

coefficients <- coef(model)
se <- summary(model)$coefficients[, "Std. Error"]
ci_lower <- coefficients - 1.96 * se
ci_upper <- coefficients + 1.96 * se

# Create a bar plot of coefficients with 95% confidence intervals

barplot(coefficients, names.arg = names(coefficients), 
        main = "Coefficients of Logistic Regression Model",
        ylab = "Coefficient Estimate", xlab = "Predictors", 
        col = ifelse(ci_lower > 0 & ci_upper > 0, "blue", 
                     ifelse(ci_lower < 0 & ci_upper < 0, "red", "gray")), 
        ylim = c(min(ci_lower), max(ci_upper)), 
        xlim = c(0.5, length(coefficients) + 0.5))

# Add error bars for confidence intervals

arrows(x0 = 1:length(coefficients), y0 = ci_lower, x1 = 1:length(coefficients), y1 = ci_upper, angle = 90, code = 3, length = 0.1)

# Add horizontal line at y = 0
abline(h = 0, lty = 2)