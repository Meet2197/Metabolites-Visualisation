library(MatchIt)
library(cobalt)
library(twang)
library(ipw)
library(caret)
library(randomForest)

metformin_psm <- ALL[, c('eid', 'Age_AC', 'Gender','BMI','Diabetes','metformin', 'MASLD','nash','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')]
metformin_psm <- subset(metformin_psm, Diabetes == 1)
pioglitazone_psm1 <- ALL_pioglitazone[, c('eid', 'Age_AC', 'Gender','BMI','Diabetes','pioglitazone', 'MASLD','nash','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')]
pioglitazone_psm <- subset(pioglitazone_psm1, Diabetes == 1)
ramipril_psm <- ALL_ramipril[, c('eid','ramipril','Age_AC', 'Gender','BMI','Diabetes', "hypertension",'MASLD','nash','alcoholicliver','toxicliver','Liverfailure','ChronHepatitis','fibrosecirrho','inflammliver','otherliverD','UnclassLiverd')]
death_psm1 <- ALL_Death [, c("eid","metformin", "Diabetes", "Age_AC", "Gender", "BMI", "liverdisease","death")]

# psm model for propensity score matching :

psm_model1 <- matchit(metformin ~ Age_AC + Gender + BMI, data = metformin_psm, method = "nearest", ratio = 5)
psm_model2 <- matchit(pioglitazone ~ Age_AC + Gender + BMI, data = pioglitazone_psm, method = "nearest", ratio = 5)
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

metformin_fs <- metformin_psm[, c("Age_AC", "Gender", "BMI")] # metformin and on metformin 
metformin_fsa <- metformin_psm$metformin

ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 10) # Perform Recursive Feature Elimination
result <- rfe(x = metformin_fs, y = metformin_fsa, sizes = c(1:3), rfeControl = ctrl)

selected_features <- predict(result$fit, metformin_fs) # Get the selected features
conf_matrix <- table(selected_features, metformin_fsa) # Calculate overall statistics for selected features
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix) # Calculate accuracy
precision <- diag(conf_matrix) / rowSums(conf_matrix) # Calculate precision
recall <- diag(conf_matrix) / colSums(conf_matrix) # Calculate recall

f1_score <- 2 * precision * recall / (precision + recall) # Calculate F1-score

cat("Confusion Matrix:\n") # Output overall statistics
print(conf_matrix)
cat("\n")
cat("Accuracy:", accuracy, "\n")
cat("Precision:", mean(precision), "\n")
cat("Recall:", mean(recall), "\n")
cat("F1-Score:", mean(f1_score), "\n")
