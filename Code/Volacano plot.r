library(PheWAS)
library(broom)
library(ggplot2)
library(ggrepel)

df <-merge(metabolites, medication, by.x="eid", all.x = TRUE )

summary_table <- function(df, metformin_var) {
  
  num_vars <- ncol(df)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 3:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    
    formula <- reformulate(c(names(x)[i], metformin_var), response = names(df)[i])
    models[[i]] <- lm(formula, data = df)
    first_tables[[i]] <- broom::tidy(models[[i]])
  } 
  
  # Combine the rows together into a final table
  final_table <- do.call("rbind", first_tables)
  
  return(final_table)
}

# Assuming 'metformin' is a column in your df data frame

final_table1 <- summary_table(df, "metformin")

final_table_Metformin <- subset(final_table1, term=="Metformin")
final_table_Metformin$logp <- -(log10(final_table_Metformin$p.value))
final_table_Metformin<-final_table_Metformin[-c(1,2,3,4,5)]  

newtable_Metformin<- aggregate(df[, 3:171], list(df$Metformin), mean, na.rm=TRUE)
newtable_Metformin[3,]=c(3,newtable_Metformin[2,-1]/newtable_Metformin[1,-1])
newtable_Metformin[4,]=log2(newtable_Metformin[3,])

newtable_Metformin<-newtable_Metformin[-c(1,2)]
table_Metformin<-as.data.frame(t(as.data.frame(newtable_Metformin)))
table_Metformin<-table_Metformin[-c(1,2,3)]
table_Metformin<-as.data.frame(table_Metformin)
all_Metformin<-cbind(table_Metformin, final_table_Metformin)

## Write the final table to your working directory as a CSV

WriteXLS::WriteXLS(all_Metformin,ExcelFileName ="Metformin_univariate.xlsx")

#Make a volcano plot 

library(ggplot2)
library(ggrepel)

# Assuming all_Metformin has 'p_value' and 'smd' columns
all_Metformin$logp <- -log10(all_Metformin$p_value) # Transform p-value
threshold <- 0.05 # Define significance threshold for p-value

ggplot(all_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = p_value < threshold)) + # Color points by significance
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + # Color significant points in red
  labs(title = "Volcano Plot of Metformin Data", x = "Standard Mean Difference (SMD)", y = "-log10(p-value)") +
  theme_minimal() # Use a minimal theme for aesthetics
