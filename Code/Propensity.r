library(MatchIt)
library(cobalt)
library(twang)
library(ipw)

subset_ALL <- ALL[, c("eid","metformin", "Diabetes", "Age_AC", "Gender", "BMI","MASLD","nash")]
subset_ALL <- subset(subset_ALL, Diabetes == 1)
pioglitazone_psm1 <- ALL_pioglitazone[, c("eid","pioglitazone", "Diabetes", "Age_AC", "Gender", "BMI","MASLD","nash")]
pioglitazone_psm <- subset(pioglitazone_psm, Diabetes == 1)
ramipril_psm <- ALL_ramipril[, c("eid","ramipril", "Diabetes", "Age_AC", "Gender", "BMI", "hypertension","MASLD","nash")]
death_psm1 <- ALL_Death [, c("eid","metformin", "Diabetes", "Age_AC", "Gender", "BMI", "liverdisease","death")]

# psm model for propensity score matching :

psm_model1 <- matchit(metformin ~ Age_AC + Gender + BMI, data = subset_ALL, method = "nearest", ratio = 5)
psm_model2 <- matchit(pioglitazone ~ Age_AC + Gender + BMI, data = pioglitazone_psm1, method = "nearest", ratio = 5)
psm_model3 <- matchit(ramipril ~  Diabetes + Age_AC + Gender + BMI + hypertension, data = ramipril_psm , method = "nearest", ratio = 5)
death_psm1 <- matchit(metformin ~ Age_AC + Gender + BMI + liverdisease + death, data = death_psm1, method = "nearest", ratio = 5)

# View summary of the matching

summary(psm_model1)
summary(psm_model2)
summary(psm_model3)
summary(death_psm1)

# Extract matched data
matched1 <- match.data(psm_model1)
matched2 <- match.data(psm_model2)
matched3 <- match.data(psm_model3)
matched4 <- match.data(psm_model3)
# trnasforming NAs

subset_ALL$BMI[is.na(subset_ALL$BMI)] <- 0
subset_ALL$BMI[!is.finite(subset_ALL$BMI)] <- 0
pioglitazone_psm1$BMI[is.na(pioglitazone_psm1$BMI)] <- 0
pioglitazone_psm1$BMI[!is.finite(pioglitazone_psm1$BMI)] <- 0
ramipril_psm$BMI[is.na(ramipril_psm$BMI)] <- 0
ramipril_psm$BMI[!is.finite(ramipril_psm$BMI)] <- 0
ramipril_psm$hypertension[!is.finite(ramipril_psm$hypertension)] <- 0
