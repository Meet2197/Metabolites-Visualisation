library(MatchIt)
library(cobalt)
library(twang)

metformin_psm <- ALL[, c('eid', 'Age_AC', 'Gender','BMI','Diabetes','metformin', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','low_fat','high_fat')]
metformin_psm <- subset(metformin_psm, Diabetes == 1)
pioglitazone_psm1 <- ALL_pioglitazone[, c('eid', 'Age_AC', 'Gender','BMI','Diabetes','pioglitazone', 'MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','low_fat','high_fat')]
pioglitazone_psm <- subset(pioglitazone_psm1, Diabetes == 1)
ramipril_psm <- ALL_ramipril[, c('eid','Age_AC', 'Gender','BMI','Diabetes','ramipril','hypertension','MASLD','MASH','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd','ALT','GLS','AST','GGT','liverdisease','low_fat','high_fat')]
death_psm1 <- ALL_Death [, c("eid","metformin", "Diabetes", "Age_AC", "Gender", "BMI", "liverdisease","death")]
normal_psm <- ALL_normal[, c('eid', 'Age_AC', 'Gender','BMI','Diabetes','metformin', 'MASLD','MASH','ALT','AST','GGT','GLS')]

# psm model for propensity score matching :

psm_model1 <- matchit(metformin ~ Age_AC + Gender + BMI + ALT + GLS + AST + GGT, data = metformin_psm, method = "nearest", ratio = 5)
psm_model2 <- matchit(pioglitazone ~ Age_AC + Gender + BMI + ALT + GLS + AST + GGT, data = pioglitazone_psm, method = "nearest", ratio = 5)
psm_model3 <- matchit(ramipril ~  Diabetes + Age_AC + Gender + BMI + hypertension + ALT + GLS + AST + GGT, data = ramipril_psm , method = "nearest", ratio = 5)
death_psm1 <- matchit(metformin ~ Age_AC + Gender + BMI + liverdisease + death, data = death_psm1, method = "nearest", ratio = 5)
psm_model5 <- matchit(metformin ~ Age_AC + Gender + BMI + ALT + GLS + AST + GGT, data = normal_psm, method = "nearest", ratio = 5)

# View summary of the matching

summary(psm_model1)
summary(psm_model2)
summary(psm_model3)
summary(death_psm1)
summary(psm_model5)

# Extract matched data
matched1 <- match.data(psm_model1)
matched2 <- match.data(psm_model2)
matched3 <- match.data(psm_model3)
matched4 <- match.data(psm_model3)
matched5 <- match.data(psm_model5)

# trnasforming NAs

metformin_psm$GGT[is.na(metformin_psm$GGT)] <- 0
metformin_psm$GGT[!is.finite(metformin_psm$GGT)] <- 0
pioglitazone_psm1$BMI[is.na(pioglitazone_psm1$BMI)] <- 0
pioglitazone_psm1$BMI[!is.finite(pioglitazone_psm1$BMI)] <- 0
ramipril_psm$BMI[is.na(ramipril_psm$BMI)] <- 0
ramipril_psm$BMI[!is.finite(ramipril_psm$BMI)] <- 0
ramipril_psm$hypertension[!is.finite(ramipril_psm$hypertension)] <- 0
