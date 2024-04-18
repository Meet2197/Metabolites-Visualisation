library(readr)
library(ggrepel)
library(EnhancedVolcano)
library(effectsize)
library(broom)
library(plotly)

# File upload : 

metabolite_metformin <-merge(matched1, metabolites, by.x="eid" , by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','metformin')

summary_table <- function(x) { # library broom is required
  num_vars <- ncol(x)
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  # Loop to create each row for the final table
  for (i in 1:num_vars) {
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_metformin)
    first_tables[[i]] <- broom::tidy(models[[i]])}  
  # Combine the rows together into a final table
  final_table <- do.call("rbind", first_tables)
  return(final_table)}
#add here all columns you want to include 
final_table1 <- summary_table(metabolite_metformin[2:169])
final_table_Metformin <- subset(final_table1, term=="metformin")

# P value generation of metformin associated metabolites against log10. 
final_table_Metformin$logp <- -(log10(final_table_Metformin$p.value))
final_table_Metformin <- final_table_Metformin[-c(1,2,3,4,5)]  

# generate smd value. 

smd_values <- apply(metabolite_metformin[, 2:169], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_metformin$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_metformin$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)}
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_df <- data.frame(smd_values)

# Adds new row to df
smd_df[3,]=c(3,smd_df[2,-1]/smd_df[1,-1])

# Add another new row at index 4
smd_df[4,] <- c(4, log2(smd_df[3,]))

# Remove the first three rows
table_smd<-as.data.frame(t(as.data.frame(smd_df)))
table_smd<-table_smd[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_Metformin<-cbind(table_smd, final_table_Metformin)

# name change for smd 
colnames(all_Metformin)[colnames(all_Metformin) == "4"] <- "smd"

smd <- all_Metformin$smd
logp <- all_Metformin$logp

# Create a hexbin plot

ggplot(all_Metformin, aes(x = smd, y = logp)) +
  geom_hex(bins = 40) + scale_fill_viridis_c() + # Optional: better color scale
  labs(title = "Hexbin Plot of metformin consuming patients", x = "SMD", y = "-log10(p-value)") + theme_minimal() +
  annotate("text", x = all_Metformin$smd, y = all_Metformin$logp, label = rownames(all_Metformin),color = "green", size = 2,
  hjust = 0.5, vjust = 0.5)

# Define significance threshold for p-value

threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 

all_Metformin <- subset(all_Metformin, logp >= 0 & logp <= 100 & smd >= -0.25 & smd <= 0.75)

# plot generation for for Metformin (this plot require ggrepl abd ggplot)

ggplot(all_Metformin, aes(x = smd, y = logp)) + geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_Metformin)),size = 3, box.padding = unit(0.25, "lines"), max.overlaps = 75) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + labs(title = "Volcano Plot of Metformin Data", x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  geom_vline(xintercept = c(-0.05, 0.05), col = "red") + geom_hline(yintercept = -log10(0.05), col = "red") + coord_cartesian(ylim = c(0, 125), xlim = c(-1 ,1)) + theme_minimal()

# fold change for enhanced volcano plot
# Set up cutoffs and parameters

ev_plot <- EnhancedVolcano(
  all_Metformin,  lab = rownames(all_Metformin),
  x = 'smd',  y = 'logp',  selectLab = rownames(all_Metformin[which(all_Metformin$logp > threshold), ]),
  pCutoff = threshold,  FCcutoff = 0.01, pointSize = 2.0, labSize = 0.0, 
  title = "Enhanced Volcano Plot of Metformin Data",subtitle = "Red points exceed threshold",
  caption = "SMD vs -log10(p-value)", xlab = "Standard Mean Difference (SMD)",
  ylab = "-log10(p-value)", xlim = c(-0.5, 1), ylim = c(-1, 1))
ev_plot + geom_text_repel(aes(label = rownames(all_Metformin)), box.padding = unit(0.35, "lines"), max.overlaps = 300)

# Create volcano plot:
volcano_plot <- plot_ly(all_Metformin, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_Metformin$smd >= -0.05 & all_Metformin$smd <= 0.05, "green", ifelse(all_Metformin$smd > 0.05, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_Metformin), "<br>", "SMD:", all_Metformin$smd, "<br>", "-log10(p-value):", all_Metformin$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_Metformin), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of Metformin Data", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
  shapes = list(list(type = "line", x0 = -0.05, x1 = -0.05, y0 = 0, y1 = 100, line = list(color = "darkblue")),
  list(type = "line", x0 = 0.05, x1 = 0.05, y0 = 0, y1 = 100, line = list(color = "darkblue")),
  list(type = "line", x0 = -0.2, x1 = 0.75, y0 = 10, y1 = 10, line = list(color = "red"))))
volcano_plot

# Diabetic Patients 

metabolite_Diabetic <-merge(matched1, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','metformin','Diabetes')

metabolite_Diabetic <- metabolite_Diabetic %>% filter(Diabetes == 1)

summary2 <- function(x) {
  num_vars <- ncol(x)
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_Diabetic)
    first_tables[[i]] <- broom::tidy(models[[i]])
    }  
  # Combine the rows together into a final table
  final_table2 <- do.call("rbind", first_tables)
  return(final_table2)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table2 <- summary2(metabolite_Diabetic[2:170])

# metformin filter with table : 

final_table_Diabetes <- subset(final_table2, term=="metformin")

# metformin filter with table : 
# P value generation of metformin associated metabolites against log10. 

final_table_Diabetes$logp <- -(log10(final_table_Diabetes$p.value))
final_table_Diabetes <- final_table_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_diabetes <- apply(metabolite_Diabetic[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_Diabetic$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_Diabetic$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds})

# df generation

smd_Diabetic_df <- data.frame(smd_diabetes)

# Adds new row to df

smd_Diabetic_df[3,]=c(3,smd_Diabetic_df[2,-1]/smd_Diabetic_df[1,-1])

# Add another new row at index 4
smd_Diabetic_df[4,] <- c(4, log2(smd_Diabetic_df[3,]))

# Remove the first three rows
table_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_Diabetic_df)))
table_smd_Diabetic<-table_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_diabetes_Metformin<-cbind(table_smd_Diabetic, final_table_Diabetes)

# name change for smd 
colnames(all_diabetes_Metformin)[colnames(all_diabetes_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 

all_diabetes_Metformin <- subset(all_diabetes_Metformin, smd >= -0.3 & smd <= 0.8 & logp >= 0 & logp <= 100)
ggplot(all_diabetes_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_diabetes_Metformin)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of Metformin consuming Patients with Diabetes",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(-20, 170)) +
  theme_minimal()

# Volcano plot: 
volcano_diabetes_Metformin <- plot_ly(all_diabetes_Metformin, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_diabetes_Metformin$smd >= -0.02 & all_diabetes_Metformin$smd <= 0.02,"green", ifelse(all_diabetes_Metformin$smd > 0.02, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_diabetes_Metformin),"<br>", "SMD:", all_diabetes_Metformin$smd, "<br>", "-log10(p-value):", all_diabetes_Metformin$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_diabetes_Metformin), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of Metformin consuming Patients with Diabetes", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.02, x1 = -0.02, y0 = 0, y1 = 100, line = list(color = "darkblue")),
           list(type = "line", x0 = 0.02, x1 = 0.02, y0 = 0, y1 = 100, line = list(color = "darkblue")),
           list(type = "line", x0 = -0.3, x1 = 0.8, y0 = 10, y1 = 10, line = list(color = "red"))))
volcano_diabetes_Metformin

# MASH Patients 

metabolite_MASLD <-merge(matched1, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','metformin','MASLD')

metabolite_MASLD <- metabolite_MASLD %>%
  filter(MASLD == 1)

summary3 <- function(x) {
  num_vars <- ncol(x)
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_MASLD)
    first_tables[[i]] <- broom::tidy(models[[i]])}  
  # Combine the rows together into a final table
  final_table3 <- do.call("rbind", first_tables)
  return(final_table3)
  final_table <- do.call("rbind", first_tables)
  return(final_table3)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 

final_table3 <- summary3(metabolite_MASLD[2:170])

# metformin filter with table : 

final_table_MASLD2 <- subset(final_table3, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_table_MASLD2$logp <- -(log10(final_table_MASLD2$p.value))
final_table_MASLD2 <- final_table_MASLD2[-c(1,2,3,4,5)]  

# metformin filter with table : 

final_table_MASLD <- subset(final_table2, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_table_MASLD$logp <- -(log10(final_table_MASLD$p.value))
final_table_MASLD <- final_table_MASLD[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD <- apply(metabolite_MASLD[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_MASLD$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_MASLD$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASLD_df <- data.frame(smd_MASLD)

# Adds new row to df

smd_MASLD_df[3,]=c(3,smd_MASLD_df[2,-1]/smd_MASLD_df[1,-1])

# Add another new row at index 4
smd_MASLD_df[4,] <- c(4, log2(smd_MASLD_df[3,]))

# Remove the first three rows
t_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_MASLD_df)))
t_smd_Diabetic<-t_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_MASLD_Metformin<-cbind(t_smd_Diabetic, final_table_MASLD)

# name change for smd 
colnames(all_MASLD_Metformin)[colnames(all_MASLD_Metformin) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 
all_MASLD_Metformin <- subset(all_MASLD_Metformin, logp >= 0 & logp <= 10 & smd >= -0.5 & smd <= 1) 
ggplot(all_MASLD_Metformin, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_MASLD_Metformin)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of Metformin consuming Patients with MASLD",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 10)) +
  theme_minimal()

# volcano plot for MASLD metformin :  
volcano_diabetes_Metformin <- plot_ly(all_MASLD_Metformin, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
marker = list(size = 7,color = ifelse(all_MASLD_Metformin$smd >= -0.05 & all_MASLD_Metformin$smd <= 0.05,"green", ifelse(all_MASLD_Metformin$smd > 0.05, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_MASLD_Metformin),"<br>", "SMD:", all_MASLD_Metformin$smd, "<br>", "-log10(p-value):", all_MASLD_Metformin$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_MASLD_Metformin), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of MASLD Patients consuming metformin", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.05, x1 = -0.05, y0 = 0, y1 = 10, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.05, x1 = 0.05, y0 = 0, y1 = 10, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.5, x1 = 1, y0 = 0.5, y1 = 0.5, line = list(color = "red"))))
volcano_diabetes_Metformin

# Diabetes patients with MASLD consuming metformin:
metabolite_MASLD_Diabetes <-merge(matched1, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','metformin','MASLD','Diabetes')
metabolite_MASLD_Diabetes <- metabolite_MASLD_Diabetes %>%
  filter(MASLD == 1, Diabetes == 1)

summary4 <- function(x) {

  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ metformin, data = metabolite_MASLD_Diabetes)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table

  final_table4 <- do.call("rbind", first_tables)
  
  return(final_table4)
}

# Assuming 'metformin' is a column in your df data frame with t test of whole dataframe: 


final_table4 <- summary4(metabolite_MASLD_Diabetes[2:170])

# metformin filter with table : 

final_MASLD_Diabetes <- subset(final_table4, term=="metformin")

# P value generation of metformin associated metabolites against log10. 

final_MASLD_Diabetes$logp <- -(log10(final_MASLD_Diabetes$p.value))
final_MASLD_Diabetes <- final_MASLD_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD_diabetes <- apply(metabolite_MASLD_Diabetes[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_MASLD_Diabetes$metformin), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_MASLD_Diabetes$metformin), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASLD_diabetes_df <- data.frame(smd_MASLD_diabetes)

# Adds new row to df

smd_MASLD_diabetes_df[3,]=c(3,smd_MASLD_diabetes_df[2,-1]/smd_MASLD_diabetes_df[1,-1])

# Add another new row at index 4
smd_MASLD_diabetes_df[4,] <- c(4, log2(smd_MASLD_diabetes_df[3,]))

# Remove the first three rows
t_smd_MASLD_Diabetic<-as.data.frame(t(as.data.frame(smd_MASLD_diabetes_df)))
t_smd_MASLD_Diabetic<-t_smd_MASLD_Diabetic[-c(1,2,3)]

# merging generated  smd and metformin tables of above code:
all_MASLD_Metformin_Diabetic<-cbind(t_smd_MASLD_Diabetic, final_MASLD_Diabetes)

# name change for smd 

colnames(all_MASLD_Metformin_Diabetic)[colnames(all_MASLD_Metformin_Diabetic) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 
all_MASLD_Metformin_Diabetic <- subset(all_MASLD_Metformin_Diabetic, logp >= 0 & logp <= 2 & smd >= -0.5 & smd <= 1) 
ggplot(all_MASLD_Metformin_Diabetic, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_MASLD_Metformin_Diabetic)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of Metformin consuming Patients with MASLD(Diabetes)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal()

# Volcano plot :
volcano_diabetes_MASLD_Metformin <- plot_ly(all_MASLD_Metformin_Diabetic, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_MASLD_Metformin_Diabetic$smd >= -0.05 & all_MASLD_Metformin_Diabetic$smd <= 0.05,"green", ifelse(all_MASLD_Metformin_Diabetic$smd > 0.05, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_MASLD_Metformin_Diabetic),"<br>", "SMD:", all_MASLD_Metformin_Diabetic$smd, "<br>", "-log10(p-value):", all_MASLD_Metformin_Diabetic$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_MASLD_Metformin_Diabetic), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 8)) %>%
  layout(title = "Volcano Plot of MASLD Patients(Diabetic) consuming metformin", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.05, x1 = -0.05, y0 = 0, y1 = 2, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.05, x1 = 0.05, y0 = 0, y1 = 2, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.5, x1 = 1, y0 = 0.1, y1 = 0.1, line = list(color = "red"))))
volcano_diabetes_MASLD_Metformin

# pioglitazone

metabolite_pioglitazone <-merge(matched2, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','pioglitazone')

# summary table loop for pioglitazone:

summary_table5 <- function(x) {

  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = metabolite_pioglitazone)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  
  final_table5 <- do.call("rbind", first_tables)
  
  return(final_table5)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 


final_table5 <- summary_table5(metabolite_pioglitazone[2:170])

# pioglitazone filter with table : 

final_table_pioglitazone <- subset(final_table5, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

final_table_pioglitazone$logp <- -(log10(final_table_pioglitazone$p.value))
final_table_pioglitazone <- final_table_pioglitazone[-c(1,2,3,4,5)]  

# generate smd value. 

smd_values_pioglitazone <- apply(metabolite_pioglitazone[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_pioglitazone$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_pioglitazone$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_pioglitazone <- data.frame(smd_values_pioglitazone)

# Adds new row to df

smd_pioglitazone[3,]=c(3,smd_pioglitazone[2,-1]/smd_pioglitazone[1,-1])

# Add another new row at index 4
smd_pioglitazone[4,] <- c(4, log2(smd_pioglitazone[3,]))

# Remove the first three rows
table_smd_pioglitazone<-as.data.frame(t(as.data.frame(smd_pioglitazone)))
table_smd_pioglitazone<-table_smd_pioglitazone[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
matched2<-cbind(table_smd_pioglitazone, final_table_pioglitazone)

# name change for smd 
colnames(matched2)[colnames(matched2) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 

matched2 <- subset(matched2, logp >= 0 & logp <= 3 & smd >= -0.2 & smd <= 0.2)
ggplot(matched2, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(matched2)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of pioglitazone consuming Patients",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal()

# Volcano plot :
volcano_pioglitazone <- plot_ly(matched2, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(matched2$smd >= -0.02 & matched2$smd <= 0.02,"green", ifelse(matched2$smd > 0.02, "blue", "red"))),
  text = ~paste("Observation:", rownames(matched2),"<br>", "SMD:", matched2$smd, "<br>", "-log10(p-value):", matched2$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(matched2), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of Patients consuming pioglitazone", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.02, x1 = -0.02, y0 = 0, y1 = 3, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.02, x1 = 0.02, y0 = 0, y1 = 3, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.2, x1 = 0.2, y0 = 0.20, y1 = 0.20, line = list(color = "red"))))
volcano_pioglitazone

# Diabetic Patients consuming pioglitazone :

Diabetic_pioglitazone <-merge(matched2, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','pioglitazone','Diabetes')

Diabetic_pioglitazone <- Diabetic_pioglitazone %>%
  filter(Diabetes == 1)

summary6 <- function(x) {

  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = Diabetic_pioglitazone)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table

  final_table6 <- do.call("rbind", first_tables)
  
  return(final_table6)

}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 


final_table6 <- summary6(Diabetic_pioglitazone[2:170])

# pioglitazone filter with table : 

pioglitazone_Diabetes <- subset(final_table6, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

pioglitazone_Diabetes$logp <- -(log10(pioglitazone_Diabetes$p.value))
pioglitazone_Diabetes <- pioglitazone_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_diabetes <- apply(Diabetic_pioglitazone[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(Diabetic_pioglitazone$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(Diabetic_pioglitazone$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_Diabetic_df <- data.frame(smd_diabetes)

# Adds new row to df

smd_Diabetic_df[3,]=c(3,smd_Diabetic_df[2,-1]/smd_Diabetic_df[1,-1])

# Add another new row at index 4
smd_Diabetic_df[4,] <- c(4, log2(smd_Diabetic_df[3,]))

# Remove the first three rows
table_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_Diabetic_df)))
table_smd_Diabetic<-table_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
all_diabetes_pioglitazone<-cbind(table_smd_Diabetic, pioglitazone_Diabetes)

# name change for smd 
colnames(all_diabetes_pioglitazone)[colnames(all_diabetes_pioglitazone) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 

all_diabetes_pioglitazone <- subset(all_diabetes_pioglitazone, logp >= 0 & logp <= 3 & smd >= -0.2 & smd <= 0.2)
ggplot(all_diabetes_pioglitazone, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_diabetes_pioglitazone)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of pioglitazone consuming Patients(Diabetic)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 5)) +
  theme_minimal()

# Volcano plot :

volcano_diabetes_Metformin <- plot_ly(all_diabetes_pioglitazone, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_diabetes_pioglitazone$smd >= -0.02 & all_diabetes_pioglitazone$smd <= 0.02,"green", ifelse(all_diabetes_pioglitazone$smd > 0.02, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_diabetes_pioglitazone),"<br>", "SMD:", all_diabetes_pioglitazone$smd, "<br>", "-log10(p-value):", all_diabetes_pioglitazone$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_diabetes_pioglitazone), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of Diabetic Patients consuming pioglitazone", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
  shapes = list(list(type = "line", x0 = -0.02, x1 = -0.02, y0 = 0, y1 = 3, line = list(color = "darkblue")),
  list(type = "line", x0 = 0.02, x1 = 0.02, y0 = 0, y1 = 3, line = list(color = "darkblue")),
  list(type = "line", x0 = -0.2, x1 = 0.25, y0 = 0.2, y1 = 0.2, line = list(color = "red"))))
volcano_diabetes_Metformin

# MASLD Patients with pioglitazone :

pioglitazone_MASLD <-merge(matched2, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','pioglitazone','MASLD')

pioglitazone_MASLD <- pioglitazone_MASLD %>%
  filter(MASLD == 1)

summary8 <- function(x) {

  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns

    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = pioglitazone_MASLD)

    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table

  final_table8 <- do.call("rbind", first_tables)
  
  return(final_table8)

}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 


final_table8 <- summary8(pioglitazone_MASLD[2:170])

# pioglitazone filter with table : 

final_table_MASLD <- subset(final_table8, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

final_table_MASLD$logp <- -(log10(final_table_MASLD$p.value))
final_table_MASLD1 <- final_table_MASLD[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD_pioglitazone <- apply(pioglitazone_MASLD[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(pioglitazone_MASLD$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(pioglitazone_MASLD$pioglitazone), sd, na.rm = TRUE)

  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_MASLD_pioglitazone <- data.frame(smd_MASLD_pioglitazone)

# Adds new row to df

smd_MASLD_pioglitazone[3,]=c(3,smd_MASLD_pioglitazone[2,-1]/smd_MASLD_pioglitazone[1,-1])

# Add another new row at index 4
smd_MASLD_pioglitazone[4,] <- c(4, log2(smd_MASLD_pioglitazone[3,]))

# Remove the first three rows
t_smd_MASLD_pioglitazone<-as.data.frame(t(as.data.frame(smd_MASLD_pioglitazone)))
t_smd_MASLD_pioglitazone<-t_smd_MASLD_pioglitazone[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
MASLD_pioglitazone<-cbind(t_smd_MASLD_pioglitazone, final_table_MASLD1)

# name change for smd 
colnames(MASLD_pioglitazone)[colnames(MASLD_pioglitazone) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 
MASLD_pioglitazone <- subset(MASLD_pioglitazone, logp >= 0 & logp <= 0.9 & smd >= -1 & smd <= 1.5)
ggplot(MASLD_pioglitazone, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(MASLD_pioglitazone)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of pioglitazone consuming Patients(MASLD)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal()

# Volcano plot
volcano_MASLD_pioglitazone <- plot_ly(MASLD_pioglitazone, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(MASLD_pioglitazone$smd >= -0.1 & MASLD_pioglitazone$smd <= 0.1,"green", ifelse(MASLD_pioglitazone$smd > 0.1, "blue", "red"))),
  text = ~paste("Observation:", rownames(MASLD_pioglitazone),"<br>", "SMD:", MASLD_pioglitazone$smd, "<br>", "-log10(p-value):", MASLD_pioglitazone$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(MASLD_pioglitazone), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of MASLD Patients consuming pioglitazone", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.1, x1 = -0.1, y0 = 0, y1 = 0.9, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.1, x1 = 0.1, y0 = 0, y1 = 0.9, line = list(color = "darkblue")),
                       list(type = "line", x0 = -1, x1 = 1.5, y0 = 0.1, y1 = 0.1, line = list(color = "red"))))
volcano_MASLD_pioglitazone

# MASLD diabetic Patients consuming pioglitazone :

metabolite_MASLD_Diabete_pioglitazone <-merge(matched2, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','pioglitazone','MASLD','Diabetes')

# MASLD and diabetes filter:

metabolite_MASLD_Diabete_pioglitazone <- metabolite_MASLD_Diabete_pioglitazone %>%
  filter(MASLD == 1, Diabetes == 1)


#Summary generation :

summary9 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns

    models[[i]] <- lm(x[[i]] ~ pioglitazone, data = metabolite_MASLD_Diabete_pioglitazone)

    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table9 <- do.call("rbind", first_tables)
  
  return(final_table9)
}

# Assuming 'pioglitazone' is a column in your df data frame with t test of whole dataframe: 


final_table9 <- summary9(metabolite_MASLD_Diabete_pioglitazone[2:170])

# pioglitazone filter with table : 

MASLD_Diabete_pioglitazone <- subset(final_table9, term=="pioglitazone")

# P value generation of pioglitazone associated metabolites against log10. 

MASLD_Diabete_pioglitazone$logp <- -(log10(MASLD_Diabete_pioglitazone$p.value))
MASLD_Diabete_pioglitazone1 <- MASLD_Diabete_pioglitazone[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASLD_diabetes <- apply(metabolite_MASLD_Diabete_pioglitazone[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_MASLD_Diabete_pioglitazone$pioglitazone), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_MASLD_Diabete_pioglitazone$pioglitazone), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASLD_diabetes_df <- data.frame(smd_MASLD_diabetes)

# Adds new row to df

smd_MASLD_diabetes_df[3,]=c(3,smd_MASLD_diabetes_df[2,-1]/smd_MASLD_diabetes_df[1,-1])

# Add another new row at index 4
smd_MASLD_diabetes_df[4,] <- c(4, log2(smd_MASLD_diabetes_df[3,]))

# Remove the first three rows
t_smd_MASLD_Diabetic<-as.data.frame(t(as.data.frame(smd_MASLD_diabetes_df)))
t_smd_MASLD_Diabetic<-t_smd_MASLD_Diabetic[-c(1,2,3)]

# merging generated  smd and pioglitazone tables of above code:
all_MASLD_pioglitazone_Diabetic<-cbind(t_smd_MASLD_Diabetic, MASLD_Diabete_pioglitazone1)

# name change for smd 
colnames(all_MASLD_pioglitazone_Diabetic)[colnames(all_MASLD_pioglitazone_Diabetic) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 
all_MASLD_pioglitazone_Diabetic <- subset(all_MASLD_pioglitazone_Diabetic, logp >= 0 & logp <= 0.90 & smd >= -0.5 & smd <= 1.5)
ggplot(all_MASLD_pioglitazone_Diabetic, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_MASLD_pioglitazone_Diabetic)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 70) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of pioglitazone consuming MASLD Patients(Diabetic)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()

# volcano plot for Pioglitazone 

volcano_MASLD_pioglitazone_diabetic <- plot_ly(all_MASLD_pioglitazone_Diabetic, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_MASLD_pioglitazone_Diabetic$smd >= -0.1 & all_MASLD_pioglitazone_Diabetic$smd <= 0.1,"green", ifelse(all_MASLD_pioglitazone_Diabetic$smd > 0.1, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_MASLD_pioglitazone_Diabetic),"<br>", "SMD:", all_MASLD_pioglitazone_Diabetic$smd, "<br>", "-log10(p-value):", all_MASLD_pioglitazone_Diabetic$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_MASLD_pioglitazone_Diabetic), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of MASLD Patients(Diabetic) consuming pioglitazone", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.1, x1 = -0.1, y0 = 0, y1 = 0.90, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.1, x1 = 0.1, y0 = 0, y1 = 0.90, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.5, x1 = 1.5, y0 = 0.1, y1 = 0.1, line = list(color = "red"))))
volcano_MASLD_pioglitazone_diabetic


# ramipril p values calculation and smd values: 

metabolite_ramipril <-merge(matched3, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','ramipril')

# summary table loop for ramipril:

summary_table10 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ ramipril, data = metabolite_ramipril)
    first_tables[[i]] <- broom::tidy(models[[i]])
    
    
  }  
  
  # Combine the rows together into a final table
  final_table10 <- do.call("rbind", first_tables)
  
  return(final_table10)
}

# Assuming 'ramipril' is a column in your df data frame with t test of whole dataframe: 

final_table10 <- summary_table10(metabolite_ramipril[2:170])

# ramipril filter with table : 

final_table_ramipril <- subset(final_table10, term=="ramipril")

# P value generation of ramipril associated metabolites against log10. 

final_table_ramipril$logp <- -(log10(final_table_ramipril$p.value))
final_table_ramipril <- final_table_ramipril[-c(1,2,3,4,5)]  

# generate smd value. 

smd_values_ramipril <- apply(metabolite_ramipril[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_ramipril$ramipril), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_ramipril$ramipril), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation
smd_ramipril <- data.frame(smd_values_ramipril)

# Adds new row to df

smd_ramipril[3,]=c(3,smd_ramipril[2,-1]/smd_ramipril[1,-1])

# Add another new row at index 4
smd_ramipril[4,] <- c(4, log2(smd_ramipril[3,]))

# Remove the first three rows
table_smd_ramipril<-as.data.frame(t(as.data.frame(smd_ramipril)))
table_smd_ramipril<-table_smd_ramipril[-c(1,2,3)]

# merging generated  smd and ramipril tables of above code:
all_ramipril<-cbind(table_smd_ramipril, final_table_ramipril)

# name change for smd 
colnames(all_ramipril)[colnames(all_ramipril) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 5

# Volacno plot generation with ggplot with log 10 p value 
all_ramipril <- subset(all_ramipril, logp >= 0 & logp <= 150 & smd >= -0.1 & smd <= 0.1)
ggplot(all_ramipril, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_ramipril)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of ramipril consuming Patients",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 350)) +
  theme_minimal()

# volcano plot for ramipril 

plot_ramipril <- plot_ly(all_ramipril, x = ~smd, y = ~logp, type = "scatter", mode = "markers", 
  marker = list(size = 7,color = ifelse(all_ramipril$smd >= -0.02 & all_ramipril$smd <= 0.02, "green", ifelse(all_ramipril$smd > 0.02, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_ramipril),"<br>", "SMD:", all_ramipril$smd, "<br>", "-log10(p-value):", all_ramipril$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_ramipril), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of Patients consuming ramipril", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.02, x1 = -0.02, y0 = 0, y1 = 150, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.02, x1 = 0.02, y0 = 0, y1 = 150, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.1, x1 = 0.1, y0 = 10, y1 = 10, line = list(color = "red"))))
plot_ramipril

# Diabetic Patients 

Diabetic_ramipril <-merge(matched3, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','ramipril','Diabetes')

Diabetic_ramipril <- Diabetic_ramipril %>%
  filter(Diabetes == 1)

summary11 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ ramipril, data = Diabetic_ramipril)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table11 <- do.call("rbind", first_tables)
  
  return(final_table11)
}

# Assuming 'ramipril' is a column in your df data frame with t test of whole dataframe: 

final_table11 <- summary11(Diabetic_ramipril[2:170])

# ramipril filter with table : 

ramipril_Diabetes <- subset(final_table11, term=="ramipril")

# P value generation of ramipril associated metabolites against log10. 

ramipril_Diabetes$logp <- -(log10(ramipril_Diabetes$p.value))
ramipril_Diabetes <- ramipril_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_diabetes <- apply(Diabetic_ramipril[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(Diabetic_ramipril$ramipril), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(Diabetic_ramipril$ramipril), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_Diabetic_df <- data.frame(smd_diabetes)

# Adds new row to df

smd_Diabetic_df[3,]=c(3,smd_Diabetic_df[2,-1]/smd_Diabetic_df[1,-1])

# Add another new row at index 4
smd_Diabetic_df[4,] <- c(4, log2(smd_Diabetic_df[3,]))

# Remove the first three rows
table_smd_Diabetic<-as.data.frame(t(as.data.frame(smd_Diabetic_df)))
table_smd_Diabetic<-table_smd_Diabetic[-c(1,2,3)]

# merging generated  smd and ramipril tables of above code:
all_diabetes_ramipril<-cbind(table_smd_Diabetic, ramipril_Diabetes)

# name change for smd 
colnames(all_diabetes_ramipril)[colnames(all_diabetes_ramipril) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.5

# Volacno plot generation with ggplot with log 10 p value 
all_diabetes_ramipril <- subset(all_diabetes_ramipril, logp >= 0 & logp <= 13 & smd >= -0.06 & smd <= 0.1)
ggplot(all_diabetes_ramipril, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_diabetes_ramipril)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of ramipril consuming Patients(Diabetes)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal()

# Plot ramipril 
plot_diabetes_ramipril <- plot_ly(all_diabetes_ramipril, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_diabetes_ramipril$smd >= -0.01 & all_diabetes_ramipril$smd <= 0.01,"green", ifelse(all_diabetes_ramipril$smd > 0.01, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_diabetes_ramipril),"<br>", "SMD:", all_diabetes_ramipril$smd, "<br>", "-log10(p-value):", all_diabetes_ramipril$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_diabetes_ramipril), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of diabetes Patients consuming ramipril ", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.01, x1 = -0.01, y0 = 0, y1 = 13, line = list(color = "darkblue")),
                       list(type = "line", x0 = 0.01, x1 = 0.01, y0 = 0, y1 = 13, line = list(color = "darkblue")),
                       list(type = "line", x0 = -0.06, x1 = 0.1, y0 = 0.5, y1 = 0.5, line = list(color = "red"))))
plot_diabetes_ramipril

# MASH Patients with ramipril :

ramipril_MASH <-merge(matched3, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','ramipril','MASH')

ramipril_MASH <- ramipril_MASH %>%
  filter(MASH == 1)

summary14 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ ramipril, data = ramipril_MASH)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table14 <- do.call("rbind", first_tables)
  
  return(final_table14)
}

# Assuming 'ramipril' is a column in your df data frame with t test of whole dataframe: 

final_table14 <- summary14(ramipril_MASH[2:170])

# ramipril filter with table : 

final_table_MASH <- subset(final_table14, term=="ramipril")

# P value generation of ramipril associated metabolites against log10. 

final_table_MASH$logp <- -(log10(final_table_MASH$p.value))
final_table_MASH <- final_table_MASH[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASH_ramipril <- apply(ramipril_MASH[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(ramipril_MASH$ramipril), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(ramipril_MASH$ramipril), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASH_ramipril <- data.frame(smd_MASH_ramipril)

# Adds new row to df

smd_MASH_ramipril[3,]=c(3,smd_MASH_ramipril[2,-1]/smd_MASH_ramipril[1,-1])

# Add another new row at index 4
smd_MASH_ramipril[4,] <- c(4, log2(smd_MASH_ramipril[3,]))

# Remove the first three rows
t_smd_MASH_ramipril<-as.data.frame(t(as.data.frame(smd_MASH_ramipril)))
t_smd_MASH_ramipril<-t_smd_MASH_ramipril[-c(1,2,3)]

# merging generated  smd and ramipril tables of above code:
MASH_ramipril<-cbind(t_smd_MASH_ramipril, final_table_MASH)

# name change for smd 
colnames(MASH_ramipril)[colnames(MASH_ramipril) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.2

# Volacno plot generation with ggplot with log 10 p value 

MASH_ramipril <- subset(MASH_ramipril, logp >= 0 & logp <= 1.7 & smd >= -0.2 & smd <= 0.6)
ggplot(MASH_ramipril, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(MASH_ramipril)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of ramipril consuming MASH Patients",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 3)) +
  theme_minimal()

# Plot ramipril 
plot_ramipril_MASH <- plot_ly(MASH_ramipril, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(MASH_ramipril$smd >= -0.05 & MASH_ramipril$smd <= 0.05,"green", ifelse(MASH_ramipril$smd > 0.05, "blue", "red"))),
  text = ~paste("Observation:", rownames(MASH_ramipril),"<br>", "SMD:", MASH_ramipril$smd, "<br>", "-log10(p-value):", MASH_ramipril$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(MASH_ramipril), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of MASH  Patients consuming ramipril ", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
  shapes = list(list(type = "line", x0 = -0.05, x1 = -0.05, y0 = 0, y1 = 1.7, line = list(color = "darkblue")),
  list(type = "line", x0 = 0.05, x1 = 0.05, y0 = 0, y1 = 1.7, line = list(color = "darkblue")),
  list(type = "line", x0 = -0.2, x1 = 0.6, y0 = 0.1, y1 = 0.1, line = list(color = "red"))))
plot_ramipril_MASH

# MASH diabetic Patients consuming ramipril :

metabolite_MASH_Diabetes <-merge(matched3, metabolites, by.x="eid",by.y="eid", all.x = TRUE, all.y = TRUE ) %>%
  select('eid','TLC','TLHC','REC','VC','CLLC','LC','HC','TLTG','TGV','TGL','TGH','TLPL','PV','PL','PH','TLEC', 'CEV','CEL','CEH','TLFC','FCV','FCL','FCH', 'TLLPL','TLLPV','TLLPH','TLL','VLDL','LDL','CHDL','ADV ','ADL','ADH','P','TGPG','TLCL',
         'PC','S','AB','A1','ABA1','TLFA','DU','O3FA','O6FA','PUFA','MUFA','SFA','LA','DHA','O3Tl','O6TL','PUTL','MUTL','STL','LATL','DCATL','PUMU','O63R','ALA','GLT','GLY','HIS','TLAA','ISO','LEU','VAL','PHA','TYR','GLS','LAC','PYR','CIT',
         'H3','ACT','ACTA','ATN','CTN','ALB','GPA','CMEL','LPEL','PCMEL','CCMEL','CECEL','FCMEL','TGEL',' VLV','TLVLV','PVLV','CVLV','CEVLV','FCVLV','TVLV','LV','TLPLV','PLV','CLV','CELV','FCLV','TGLV','MV','TLPMV','PMV','CMV','CEMV','FCMV',
         'TMV','SV','TLPSV','PSV','CSV','CESV','FCSV','TGSV','VSV','TLPVSV','PVSV','CVSV','CEVSV','FCVSV','TGVSV','I','TLPI','PI','CI','CEI','FCI','TGI','LL','TLPLL','PLL','CLL','CELL','FCLL','TGLL','ML','TLPML','PML','CML','CEML','FCML',
         'TGML','SL','TLPSL','PSL','CSL','CESL','FCSL','TGSL','VLH','TLPVLH','PVLH','CVLH','CEVLH','FCVLH','TGVLH','LH','TLPLH','PLH','CLH','CELH','FCLH','TGLH','MH','TLPMH','PMH','CMHD','ramipril','MASH','Diabetes')

# MASH and diabetes filter:

metabolite_MASH_Diabetes <- metabolite_MASH_Diabetes %>%
  filter(MASH == 1, Diabetes == 1)

#Summary generation :

summary15 <- function(x) {
  
  num_vars <- ncol(x)
  
  # Pre-define lists that will be populated and then collapsed by rest of function
  models <- vector("list", length = num_vars)
  first_tables <- vector("list", length = num_vars)
  second_tables <- vector("list", length = num_vars)
  
  # Loop to create each row for the final table
  for (i in 1:num_vars) {  # Start from column 3, as columns 1 and 2 are assumed to be non-predictor columns
    models[[i]] <- lm(x[[i]] ~ ramipril, data = metabolite_MASH_Diabetes)
    first_tables[[i]] <- broom::tidy(models[[i]])
  }  
  
  # Combine the rows together into a final table
  final_table15 <- do.call("rbind", first_tables)
  
  return(final_table15)
}

# Assuming 'ramipril' is a column in your df data frame with t test of whole dataframe: 

final_table15 <- summary15(metabolite_MASH_Diabetes[2:170])

# ramipril filter with table : 

final_MASH_Diabetes <- subset(final_table15, term=="ramipril")

# P value generation of pioglitazone associated metabolites against log10. 


final_MASH_Diabetes$logp <- -(log10(final_MASH_Diabetes$p.value))
final_MASH_Diabetes <- final_MASH_Diabetes[-c(1,2,3,4,5)]  

# generate smd value. 

smd_MASH_ramipril <- apply(metabolite_MASH_Diabetes[, 2:170], 2, function(x) {
  group_means <- aggregate(x, list(metabolite_MASH_Diabetes$ramipril), mean, na.rm = TRUE)
  group_sds <- aggregate(x, list(metabolite_MASH_Diabetes$ramipril), sd, na.rm = TRUE)
  
  # Extract standard deviations excluding the first column (group indicator)
  group_sds_values <- group_sds[, -1]
  
  # Square the standard deviations
  group_sds_squared <- group_sds_values^2
  
  # Ensure proper formatting
  if (length(dim(group_sds_squared)) < 2) {
    group_sds_squared <- matrix(group_sds_squared, nrow = length(group_sds_squared), ncol = 1)
  }
  
  # Calculate row sums of squared standard deviations
  row_sums_sds_squared <- rowSums(group_sds_squared)
  
  # Divide by 2 as required for SMD calculation
  group_smds <- abs(diff(group_means$x)) / sqrt(row_sums_sds_squared / 2)
  group_smds
})

# df generation

smd_MASH_diabetes_df <- data.frame(smd_MASH_ramipril)

# Adds new row to df

smd_MASH_diabetes_df[3,]=c(3,smd_MASH_diabetes_df[2,-1]/smd_MASH_diabetes_df[1,-1])

# Add another new row at index 4
smd_MASH_diabetes_df[4,] <- c(4, log2(smd_MASH_diabetes_df[3,]))

# Remove the first three rows
t_smd_MASH_Diabetic<-as.data.frame(t(as.data.frame(smd_MASH_diabetes_df)))
t_smd_MASH_Diabetic<-t_smd_MASH_Diabetic[-c(1,2,3)]


# merging generated  smd and ramipril tables of above code:
all_MASH_ramipril_Diabetic<-cbind(t_smd_MASH_Diabetic, final_MASH_Diabetes)

# name change for smd 
colnames(all_MASH_ramipril_Diabetic)[colnames(all_MASH_ramipril_Diabetic) == "4"] <- "smd"

# Define significance threshold for p-value
threshold <- 0.05

# Volacno plot generation with ggplot with log 10 p value 

all_MASH_ramipril_Diabetic <- subset(all_MASH_ramipril_Diabetic, logp >= 0 & logp <= 2 & smd  >= -0.3 & smd <= 1.5)
ggplot(all_MASH_ramipril_Diabetic, aes(x = smd, y = logp)) +
  geom_point(aes(color = logp > threshold)) +
  geom_text_repel(aes(label = rownames(all_MASH_ramipril_Diabetic)),size = 2, box.padding = unit(0.25, "lines"), max.overlaps = 35) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Volcano Plot of ramipril consuming MASH Patients(Diabetic)",x = "Standard Mean Difference (SMD)",y = "-log10(p-value)") +
  coord_cartesian(ylim = c(0, 3)) +
  theme_minimal()

# Plot Ramipril : 

plot_diabetes_ramipril_MASH2 <- plot_ly(all_MASH_ramipril_Diabetic, x = ~smd, y = ~logp, type = "scatter", mode = "markers",
  marker = list(size = 7,color = ifelse(all_MASH_ramipril_Diabetic$smd >= -0.05 & all_MASH_ramipril_Diabetic$smd <= 0.05,"green", ifelse(all_MASH_ramipril_Diabetic$smd > 0.05, "blue", "red"))),
  text = ~paste("Observation:", rownames(all_MASH_ramipril_Diabetic),"<br>", "SMD:", all_MASH_ramipril_Diabetic$smd, "<br>", "-log10(p-value):", all_MASH_ramipril_Diabetic$logp),showlegend = FALSE) %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers",  marker = list(size = 10, color = "green"), name = "SMD between -0.05 and 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "blue"), name = "SMD > 0.05") %>%
  add_trace(x = NA, y = NA, type = "scatter", mode = "markers", marker = list(size = 10, color = "red"), name = "SMD < -0.05") %>%
  add_trace(type = "scatter", mode = "text", text = ~rownames(all_MASH_ramipril_Diabetic), hoverinfo = "text", textposition = "top center", showlegend = FALSE, textfont = list(size = 10)) %>%
  layout(title = "Volcano Plot of MASH (Diabetic) Patients consuming ramipril ", xaxis = list(title = "Standard Mean Difference (SMD)"), yaxis = list(title = "-log10(p-value)"),
         shapes = list(list(type = "line", x0 = -0.05, x1 = -0.05, y0 = 0, y1 = 2, line = list(color = "darkblue")),
         list(type = "line", x0 = 0.05, x1 = 0.05, y0 = 0, y1 = 2, line = list(color = "darkblue")),
         list(type = "line", x0 = -0.3, x1 = 1.5, y0 = 0.1, y1 = 0.1, line = list(color = "red"))))
plot_diabetes_ramipril_MASH2
