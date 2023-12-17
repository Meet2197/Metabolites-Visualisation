# File for feature selection Package:

library(glmnet)

Feature_selection <- read.csv(file = "C:/Users/User/Desktop/PhD Documentation/My drafts/Feature_selection.csv")

# File rename for gene selection:

setnames(Feature_selection, old=colnames(Feature_selection), new = c('Eid','sex','BMI','diabetes','HbA1c','ALT','AST','GGT','T_C','HDL_C','LDL_C','SBP','DBP','smoking','waist_cir','hip_cir','Phy_act','frt_int','veg_int','fish_int','meat_int','hist_liver','trg','fsb','ethb','edu','occu','income','Apo_A','Apo_B','Uric','Meno','Hormon_Re','Ultra','Mri','ICD_code'))

# Feature_selection<- fread("C:/Users/User/Desktop/Data/ukb52200.csv", select=c("eid","31-0.0","21001-0.0","2443-0.0","30750-0.0","30690-0.0","30680-0.0","30710-0.0","30760-0.0","30780-0.0","30740-0.0","4080-0.0","4079-0.0","20116-0.0","23104-0.0","23105-0.0","864-0.0","1031-0.0","103990-0.0","103140-0.0","1349-0.0","20107-0.0","30810-0.0","30770-0.0","21000-0.0","6138-0.0","22601-0.0","738-0.0","30600-0.0","30610-0.0","30790-0.0","2724-0.0","20003-0.0","22400-0.0","21022-0.0","41270-0.0"))
# write.csv(Feature_selection, "C:/Users/User/Desktop/PhD Documentation/My drafts/Feature_selection.csv", row.names=FALSE)

# Data training for feature selection: 

sample_index <- sample(nrow(Feature_selection), 0.8 * nrow(Feature_selection))
train_data <- Feature_selection[sample_index, ]
test_data <- Feature_selection[-sample_index, ]

# Matrix data for training data: 

x_train <- as.matrix(train_data[, -1])  # Exclude the 'Eid' column
y_train <- train_data$  # Replace 'YourOutcomeVariable' with your actual outcome variable name

# Complete cases for training data:    
  
complete_cases <- complete.cases(x_train, y_train)
x_train <- x_train[complete_cases, ]
y_train <- y_train[complete_cases]
  
# Lasso Model 

lasso_model <- glmnet(x_train, y_train, alpha = 1) 

response_var <- Feature_selection$diabetes
predictors <- Feature_selection[, -c(1)]  # Exclude 'Eid' and 'BMI'

if (nrow(predictors) != length(response_var)) {
  stop("Number of rows in 'predictors' does not match the length of 'response_var'.")
}

# Handle missing values, if any 

predictors <- na.omit(predictors)
response_var <- na.omit(response_var)

lasso_model <- cv.glmnet(as.matrix(predictors), response_var, alpha=1)

# Plot the cross-validation results to find the optimal lambda

plot(lasso_model)