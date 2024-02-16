library(MatchIt)
library(cobalt)


subset_ALL <- ALL[, c("metformin", "Diabetes", "Age_AC", "Gender", "BMI")]
subset_ALL_diabetes <- subset(subset_ALL, Diabetes == 1)

# psm model for propensity score matching :

psm_model <- matchit(metformin ~ Diabetes + Age_AC + Gender + BMI, data = ALL, method = "nearest", ratio = 5)

# View summary of the matching
summary(psm_model)

# Extract matched data
matched_data <- match.data(psm_model)
psm_model_df <- data.frame(matched_data)