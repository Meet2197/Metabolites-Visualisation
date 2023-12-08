#Libraries used in these file 

library(plotly)
library(DescTools)
library(dplyr)
library(knitr)
library(tidyverse)
library(naniar)
library(data.table)
library(readr)
library(lubridate)
library(viridis)
library(DT)
library(arsenal)
library(car)

#File formation here as d,e,f,g are new division of file to run as data set to read csv format file 

d = read.csv("C:\\Users\\Bhatt\\Desktop\\M2_internsip\\data_files\\Aliquots.csv", sep = ";")
e = read.csv("C:\\Users\\Bhatt\\Desktop\\M2_internsip\\data_files\\Patients.csv", sep = ";")
f = read.csv("C:\\Users\\Bhatt\\Desktop\\M2_internsip\\data_files\\Samples.csv", sep = ";")
g = read.csv("C:\\Users\\Bhatt\\Desktop\\M2_internsip\\data_files\\Visits.csv", sep = ";")

#Null data cleaning from available data set d[d == "NULL"] <- NA, e[e == "NULL"] <- NA , f[f == "NULL"] <- NA ,g[g == "NULL"] <- NA

#new column names as per data dictionary

setnames(d, old=("DZIF_TEMPERATURE"), new = c('Temperature'))
setnames(e, old=colnames(e), new = c('ID','Code','Study','Gender','Birthyear','Affilation TTU','Loss-to-follow-up','Informed Consent','Consent Restriction','Withdrawl'))
setnames(f, old=colnames(f), new = c('SAMPLEID','CXXSample-ID','PatientId','VisitID','SAMPLETYPE','SAMPLINGDATE','INITIALAMOUNT','INITIALUNIT','AMOUNTREST','RESTUNIT','REPOSITIONDATE','SAMPLELOCATION','Processing  Time','Sampling location','OTHER SAMPLE TYPE','Alkohol','Fasting status','anesthesia','Deviations blood collection','Anuria','Deviations urine','Menstrual bleeding','Deviations from extraction','Storage temperature sample','Longterm storage sample'))
setnames(g, old=colnames(g), new = c('PATIENTID','EPISODEID','Sampling Time','Processing Time','Visit-Type','VISIT DETAILS','Diagnoses','DIAGDATE','COMMENTS','Secondary DIAGNOSISCODE','Secondary DIAGDATE',' Secondary COMMENTS','Height','Weight','Smoking status','New Infection','Type of infection','Prescription antibiotics','Outpatient medical treatment','Inpatient medical treatment','SURGERY'))

#Values change

# myValues = list(list("BOR_TBC", "Borstel Tuberkulose"), list("H_TX", "Hannover"), list("HD_TX", "Heidelberg"), list("HHBNI_FWS", "Hamburg FWS"), list("HHBNI_FWSC", "HamburgFWS children"), list("K_HIV", "KOLN HIV"), list("MLMU_TX","Munchen(LMU)"),list("MTUM_TX","Munchen (TUM) TX"),list("TUE_TX","Tubingen TX"))for (values in myValues){e$Study[e$Study == as.character(values[1])] = as.character(values[2])}
myValues = list(list("Null", "Unlocated"))
for (values in myValues){
  e$Study[e$Study == as.character(values[1])] = as.character(values[2])
}
for (values in myValues){
  d$SAMPLETYPE[d$SAMPLETYPE == as.character(values[1])] = as.character(values[2])
}
myValues = list(list("LIQUID_CSF", "CSF"), list("LIQUID_EDTA", "EDTA"), list("LIQUID_STO_NAT", "Stuhl-nativ"), list("LIQUID_STO_STA", "Stuhl-stabilisiert"), list("LIQUID_SER", "SERUM"), list("LIQUID_RNA_PAX", "Paxgene RNA"), list("LIQUID_URI","URINE"),list("LIQUID_DNA","DNA"),list("LIQUID_BC","Buffy Coat"),list("LIQUID_RNA","RNA"),list("LIQUID_PLA","Plasma"),list("LIQUID_URI_SED","Urin-Sediment"),list("TISSUE_PATH","Pathogen Tissue"),list("LIQUID_PBMC","PBMCs")) 

for (values in myValues){
  d$SAMPLETYPE[d$SAMPLETYPE == as.character(values[1])] = as.character(values[2])
}
myValues = list(list("LIQUID_CSFC", "CSFC"), list("LIQUID_WB", "WB"), list("LIQUID_STO_NAT", "Stuhl-nativ"), list("LIQUID_STO_STA", "Stuhl-stabilisiert"), list("LIQUID_SER", "SERUM"), list("LIQUID_RNA_PAX", "Paxgene RNA"), list("LIQUID_URI","URINE"),list("LIQUID_DNA","DNA"),list("LIQUID_BC","Buffy Coat"),list("LIQUID_RNA","RNA"),list("LIQUID_EDTA","EDTA"),list("UBERKANNT","Undefined"),list("LIQUID_CSF","CSF"),list("LIQUID_OTH","OTH"),list("LIQUID_TS","TS")) 

for (values in myValues){
  f$SAMPLETYPE[f$SAMPLETYPE == as.character(values[1])] = as.character(values[2])
}
myValues = list(list("Default user workplace (admin)", "Admin"), list("H", "Hannover"), list("MHZM", "Munchen"), list("TUE", "Tubingen"), list("NULL","Un-categorized"))
for (values in myValues){
  f$SAMPLELOCATION[f$SAMPLELOCATION == as.character(values[1])] = as.character(values[2])
}
myValues = list(list("MAL","Malaria"), list("APS", "African Partner Sites"), list("HEP", "Hepatitis"), list("TBC","Tuberkulose"))
for (values in myValues){
  e$`Affilation TTU`[e$`Affilation TTU` == as.character(values[1])] = as.character(values[2])
}
myValues = list(list("H","Hannover"), list("HD","Heidelberg"), list("HHBNI", "Hamburg (BNI)"), list("K", "Koln"), list("TUE","Tubingen"), list("NULL","Un-categorized"))
for (values in myValues){
  d$SAMPLELOCATION[d$SAMPLELOCATION == as.character(values[1])] = as.character(values[2])
}

#Removing Charachter from each row in one column 

g$`VISIT DETAILS` = gsub("VisitID: ","",g$`VISIT DETAILS`)
g$`VISIT DETAILS` = gsub(", Follow-up 28 days after hospital admission","",g$`VISIT DETAILS`)
g$'Visit-Type' = gsub("DZIF_","",g$'Visit-Type')

# Numeric format change

g$`VISIT DETAILS` <- as.numeric(g$`VISIT DETAILS`)
e$Code <- as.numeric(e$Code)
e$`Affilation TTU` = gsub("TTU_","",e$`Affilation TTU`)
g$PATIENTID <- as.numeric (g$PATIENTID)
g$EPISODEID <- as.numeric(g$EPISODEID)
e$Code <- as.numeric(e$Code)
g$`VISIT DETAILS` <- as.numeric(g$`VISIT DETAILS`)
e$ID <- as.numeric(e$ID)
f$SAMPLEID <- as.numeric(f$SAMPLEID)
as.numeric(g$Height)
as.numeric(g$Weight)

# Split in one column values:

for (code in e$Code){
  x = unlist(str_split(code, "_"))
  e$Code[e$Code == code] = x[length(x)]}

# Freezing time 

g$`Sampling Time` <- as.numeric(g$`Sampling Time`)
g$`Processing Time` <- as.numeric(g$`Processing Time`)
g$`Freezing duration` <- with(g,(difftime(`Sampling Time`, `Processing Time`) ))

# Combining two columns 

common <- intersect(e$ID, g$PATIENTID)  
e_common = e[common,] # give you common rows in data frame 1  
g_common = g[common,]
common_comb = merge(e, g, by.x ="ID" , by.y ="PATIENTID") %>%
  select(Study, `Affilation TTU`,Height ,Weight ,'Smoking status','Sampling Time','Processing Time','Visit-Type', Diagnoses)
for (col in unique(common_comb$Study)){
  print(common_comb[common_comb$Study == col ,])}

## Missing values of studies ##

n_study = length(unique(common_comb$Study))
hgroupDF <- as.data.frame(colnames(common_comb))
names(hgroupDF) <- c("Parameters")
for(s in unique(common_comb$Study)){
  h = common_comb[which(common_comb$Study == s),]
  h[h == "NULL"] = NA
  h_count = sum(is.na(h))
  columnList = colnames(h)
  h_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
  hcountDF <- as.data.frame(matrix(ncol=2, nrow=0))
  for(c in 0:length(columnList)){
    hcountDF[c,1] = as.character(columnList[c])
    hcountDF[c,2] = (sum(is.na(h[columnList[c]]))/count(h))*100 
  }
  hgroupDF[[s]] <- hcountDF$V2
}
hcountDF <- hgroupDF %>% gather(colname, value, -Parameters)
names(hcountDF) <- c("Parameters", "Study", "NullCount")
hcountDF %>% ggplot(mapping = aes(y = Parameters, x = Study))+geom_tile(aes(fill = Study)) +  geom_text(aes(label = format(round(NullCount , 2), nsmall = 2), ))  +ggtitle("Missing Values of Studies") +  theme(plot.title = element_text(hjust = 0.5))  + scale_fill_viridis_d(option = "D") 
hcountDF %>%   ggplot(aes(Parameters,Study, fill = Study, label = Study)) +geom_tile() + facet_wrap(~Value, scales = "Parameters") +  geom_label(color = "black", fill = "white") +scale_fill_viridis_c(guide = "none") +labs(x = "Study", y = "Parameters",title = paste0("Missing Values of Studies", titel, "Missing values")) +theme(strip.background = element_rect(fill = "#001540"),strip.text = element_text(colour = 'white', face = "bold"), legend.position = "bottom")
  
## Missing values of sites ##

hagroupDF <- as.data.frame(colnames(common_comb))
names(hagroupDF) <- c("Parameters")
for(s in unique(common_comb$`Affilation TTU`)){
  ha = common_comb[which(common_comb$`Affilation TTU` == s),]
  ha[ha == "NULL"] = NA
  ha_count = sum(is.na(ha))
  columnList = colnames(ha)
  ha_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
  hacountDF <- as.data.frame(matrix(ncol=2, nrow=0))
  for(c in 0:length(columnList)){
    hacountDF[c,1] = as.character(columnList[c])
    hacountDF[c,2] = (sum(is.na(ha[columnList[c]]))/count(ha))*100 ## sum(is.na(h[columnList[c]]))
  }
  hagroupDF[[s]] <- hacountDF$V2
}
hacountDF <- hagroupDF %>% gather(colname, value, -Parameters)
names(hacountDF) <- c("Parameters", "Sites", "NullCount")
hacountDF %>% ggplot(mapping = aes(y = Parameters, x = Sites)) +geom_tile(aes(fill = Sites)) +  geom_text(aes(label = format(round(NullCount , 2), nsmall = 2))) +ggtitle("Missing Values at partner sites") +  theme(plot.title = element_text(hjust = 0.5))  + scale_fill_viridis_d(option = "D") 
plot_ly(type = 'table',columnwidth = c(50, 50,50),columnorder = c(0, 1,2), header = list(values = c("Parameters","Sites","Missing values"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(hacountDF$Parameters,hacountDF$Sites,hacountDF$`Missing values`), align = c("center","center","center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))

# Humburg FWS

h = common_comb[which(common_comb$Study == "Hamburg FWS"),]
h[h == "NULL"] = NA
sum(is.na(h))
nulls = is.na(h)
count(h)*ncol(h)
h_count = sum(is.na(h))
h_percent <- h_count / (count(h)*ncol(h))
h1_percent <- format(round(h_percent , 2), nsmall = 2)
sum(is.na(h))/count(h)
columnList = colnames(h)[0:ncol(h)]
h_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
hcountDF <- as.data.frame(matrix(ncol=3, nrow=0))
for(c in 0:length(columnList)){
  hcountDF[c,1] = as.character(columnList[c])
  hcountDF[c,2] = sum(is.na(h[columnList[c]]))
  hcountDF[c,3] = (sum(is.na(h[columnList[c]]))/count(h))*100
}
names(hcountDF) <- c("Standard", "Null count","TTU_MAL")
hcountDF$Standard = gsub("Study","Hamburg FWS",hcountDF$Standard)
plot_ly(type = 'table',columnwidth = c(50, 50,50),columnorder = c(0, 1,2 ), header = list(values = c("Parameters","Missing values","Percentage (%)"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(hcountDF$Standard,hcountDF$`Null count`,hcountDF$percentage), align = c("center", "center", "center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))
hcountDF %>% ggplot(mapping = aes(x = TTU_MAL, y = Standard)) +geom_tile(aes(fill = Standard)) +  geom_text(aes(label = format(round(`TTU_MAL` , 2), nsmall = 2))) +ggtitle("Missing Values at Humburg FWS") +  theme(plot.title = element_text(hjust = 0.5))  + scale_fill_viridis_d(option = "D") 

#expected heatmap

   heatmap_probe <- df_probenart %>% 
      ggplot(aes(value, Probenart, fill = n, label = n)) +
      geom_tile() +
      facet_wrap(~key, scales = "free_x") +
      geom_label(color = "black", fill = "white") +
      scale_fill_viridis_c(guide = "none") +
      labs(x = "", 
           y = "",
           title = paste0("Absolute Verteilung der ", titel, "proben")) +
      theme(strip.background = element_rect(fill = "#001540"),
            strip.text = element_text(colour = 'white', face = "bold"),
            legend.position = "bottom")
 



#Humburg FWS children

hew = common_comb[common_comb$Study]
hew = common_comb[which(common_comb$Study == "Hamburg FWS children"),]
hew[,names(hew):=lapply(.SD, recode, '"NULL"=NA')]
hew[hew == "NULL"] = NA
sum(is.na(hew))
nulls = is.na(hew)
count(hew)*ncol(hew)
hew_count = sum(is.na(hew))
hew_percent <- hew_count / (count(hew)*ncol(hew))
hew_percent  
sum(is.na(hew))/count(hew)
columnList = colnames(hew)[0:ncol(hew)]
hew_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
hew_count <- hew(hew_count)
hewcountDF$Standard = gsub("Study","Hamburg FWS children",hewcountDF$Standard)
hewcountDF <- as.data.frame(matrix(ncol=3, nrow=0))
for(c in 0:length(columnList)) {
  hewcountDF[c,1] = as.character(columnList[c])
  hewcountDF[c,2] = sum(is.na(hew[columnList[c]]))
  hewcountDF[c,3] = (sum(is.na(hew[columnList[c]]))/count(hew))*100
}
names(hewcountDF) <- c("Standard", "Null count","TTU_MAL")
hewcountDF %>% ggplot(mapping = aes(x = TTU_MAL, y = Standard)) +geom_tile(mapping = aes(fill = Standard))+geom_text(aes(label = format(round(`TTU_MAL` , 2), nsmall = 2)))  +ggtitle("Missing Values of Humburg FWS Children") + theme_light() +theme(plot.title = element_text(hjust = 0.5)) + scale_fill_viridis_d(option = "D")
plot_ly(type = 'table',columnwidth = c(50, 25, 50),columnorder = c(0, 1 , 2), header = list(values = c("Parameters","Missing values","Percentage (%)"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(hewcountDF$Standard,hewcountDF$`Null count`,hewcountDF$percentage), align = c("center", "center", "center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))
hew$`Affilation TTU`

# Munchen(TUM)TX

TUM = common_comb[common_comb$Study]
TUM = common_comb[which(common_comb$Study == "MTUM"),]
TUM[,names(TUM):=lapply(.SD, recode, '"NULL"=NA')]
TUM[TUM == "NULL"] = NA
TUM_count = sum(is.na(TUM))
TUM_percent <- TUM_count / (count(TUM)*ncol(TUM))
TUM_percent  
sum(is.na(TUM))/count(TUM)
columnList = colnames(TUM)[0:ncol(TUM)]

TUM_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
TUM_count <- TUM(TUM_count)
TUMcountDF$Standard = gsub("Study","MTUM",TUMcountDF$Standard)
TUMcountDF <- as.data.frame(matrix(ncol=3, nrow=0))
for(c in 0:length(columnList)) {
  
  TUMcountDF[c,1] = as.character(columnList[c])
  TUMcountDF[c,2] = sum(is.na(TUM[columnList[c]]))
  TUMcountDF[c,3] = (sum(is.na(TUM[columnList[c]]))/count(TUM))*100
}
names(TUMcountDF) <- c("Standard", "Null count","TTU_IICH")
plot_ly(type = 'table',columnwidth = c(50, 25, 50),columnorder = c(0, 1 , 2), header = list(values = c("Parameters","Missing Values","Percentage (%)"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(TUMcountDF$Standard,TUMcountDF$`Null count`,TUMcountDF$percentage), align = c("center", "center", "center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))
TUMcountDF %>% ggplot(mapping = aes(x = TTU_IICH, y = Standard)) +geom_tile(mapping = aes(fill = Standard)) +geom_text(aes(label = format(round(`TTU_IICH` , 2), nsmall = 2))) +ggtitle("Missing values Munchen (TUM) TX") + theme_light() +theme(plot.title = element_text(hjust = 0.5)) + scale_fill_viridis_d(option = "D") 
TUM$`Affilation TTU`

# Munchen(LMU)TX

LMU = common_comb[common_comb$Study]
LMU = common_comb[which(common_comb$Study == "Munchen(LMU)TX"),]

LMU[,names(LMU):=lapply(.SD, recode, '"NULL"=NA')]

LMU[LMU == "NULL"] = NA

LMU_count = sum(is.na(LMU))
LMU_percent <- TUM_count / (count(LMU)*ncol(LMU))
LMU_percent  
sum(is.na(LMU))/count(LMU)
columnList = colnames(LMU)[0:ncol(LMU)]

LMU_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
LMU_count <- LMU(LMU_count)
LMUcountDF$Standard = gsub("Study","Munchen(LMU)TX",LMUcountDF$Standard)
LMUcountDF <- as.data.frame(matrix(ncol=3, nrow=0))
for(c in 0:length(columnList)) {
  
  LMUcountDF[c,1] = as.character(columnList[c])
  LMUcountDF[c,2] = sum(is.na(LMU[columnList[c]]))
  LMUcountDF[c,3] = (sum(is.na(LMU[columnList[c]]))/count(LMU))*100
}
names(LMUcountDF) <- c("Standard", "Null count","percentage")
plot_ly(type = 'table',columnwidth = c(50, 25, 50),columnorder = c(0, 1 , 2), header = list(values = c("Parameters","Null Count","Percentage (%)"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(TUMcountDF$Standard,TUMcountDF$`Null count`,TUMcountDF$percentage), align = c("center", "center", "center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))

LMUcountDF %>% 
  ggplot(mapping = aes(x = percentage, y = Standard)) +geom_tile(mapping = aes(fill = Standard))+geom_text(aes(label = `Null count`)) +ggtitle("Missing values of Munchen (LMU) TX") + theme_light() +theme(plot.title = element_text(hjust = 0.5))  + scale_fill_viridis_d(option = "D")

#Koln HIV

KH = common_comb[common_comb$Study]
KH = common_comb[which(common_comb$Study == "KOLN HIV"),]
KH[,names(KH):=lapply(.SD, recode, '"NULL"=NA')]
KH[KH == "NULL"] = NA
KH_count = sum(is.na(KH))
KH_percent <- KH_count / (count(KH)*ncol(KH))
KH_percent  
sum(is.na(KH))/count(KH)
columnList = colnames(KH)[0:ncol(KH)]

KH_count <-sapply(x, function(y) sum(length(which(is.na(y)))))
KH_count <- KH(KH_count)
KHcountDF$Standard = gsub("Study","KOLN HIV",KHcountDF$Standard)
KHcountDF <- as.data.frame(matrix(ncol=3, nrow=0))
for(c in 0:length(columnList)) {
  KHcountDF[c,1] = as.character(columnList[c])
  KHcountDF[c,2] = sum(is.na(KH[columnList[c]]))
  KHcountDF[c,3] = (sum(is.na(KH[columnList[c]]))/count(KH))*100
}
names(KHcountDF) <- c("Standard", "Null count","KolnHiv")
plot_ly(type = 'table',columnwidth = c(50, 25, 50),columnorder = c(0, 1 , 2), header = list(values = c("Parameters","Missing values","KolnHiv (%)"),line = list(color = '#506784'), fill = list(color = '#119DFF'), align = c('center','center','center'), font = list(color = 'white', size = 12)),cells = list(values = rbind(KHcountDF$Standard,KHcountDF$`Null count`,KHcountDF$percentage), align = c("center", "center", "center"), line = list(color = '#506784'),fill = list(color = c('#25FEFD', 'white')),font = list(family = "Arial", size = 12, color = c('#506784'))))
KHcountDF %>% ggplot(mapping = aes(x = KolnHiv, y = Standard)) +geom_tile(mapping = aes(fill = Standard))+geom_text(aes(label = format(round(KolnHiv , 2), nsmall = 2))) +ggtitle("Missing values for KOLN HIV study") + theme_light() +theme(plot.title = element_text(hjust = 0.5))  + scale_fill_viridis_d(option = "D")
KH$`Affilation TTU`

#Plot generation
DataTable(d$SAMPLETYPE,d$SAMPLELOCATION)
table(d$SAMPLETYPE,d$SAMPLELOCATION)

ALIQUOTS = table(d$SAMPLETYPE,d$SAMPLELOCATION)
COHORT = table(e$Study,e$`Affilation TTU`)
SAMPLE = table(f$SAMPLETYPE,f$SAMPLELOCATION)
STUDY = table(e$Study,e$`Affilation TTU`)

ggplot(e %>% count(Study, `Affilation TTU`), aes(fill=Study, x=n, y=`Affilation TTU`)) + geom_bar(position="stack", stat="identity") +ggtitle("Studies against Partner sites") + geom_col() + scale_fill_viridis_d(option = "D")
ggplot(f %>% count(SAMPLETYPE, SAMPLELOCATION), aes(fill=SAMPLETYPE, x=n, y=SAMPLELOCATION))  + geom_bar(position="stack", stat="identity") +ggtitle("Different samples at Partner sites") +  theme(plot.title = element_text(hjust = 0.5)) + geom_col() + scale_fill_viridis_d(option = "D")
ggplot(d %>% count(SAMPLETYPE, SAMPLELOCATION), aes(fill=SAMPLETYPE, x=n, y=SAMPLELOCATION))  + geom_bar(position="stack", stat="identity") +ggtitle("Different samples at Partner sites") +  theme(plot.title = element_text(hjust = 0.5)) + geom_col() + scale_fill_viridis_d(option = "D")

plot_ly(type = 'table', header = list(values = c("<b>Aliquots</b>","<b>Hamburg (BNI)","</b><b>Hannover</b>","<b>Heidelberg</b>","<b>Koln</b>","<b>MHZM</b>","<b>NULL</b>","<b>Tubingen</b>",names(ALIQUOTS)),line = list(color = '#506784'),  fill = list(color = 'rgb(235, 100,100,100,100, 230)'), align = c('Center', rep('center', ncol(ALIQUOTS))), font = list(color = 'white', size = 12)),cells = list(values = rbind(rownames(ALIQUOTS),t(as.matrix(unname(ALIQUOTS)))), line = list(color = '#506784'),fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222,222, 249, 0.65)')),font = list(family = "Arial", size = 12, color = c('#506784'))))
plot_ly(type = 'table', header = list(values = c("<b>Studies</b>","<b>IICH","</b><b>MAL</b>","<b>NULL</b>","<b>TBC</b>",names(SAMPLE)),line = list(color = '#506784'),  fill = list(color = 'rgb(235, 100,100,100,100, 230)'), align = c('Center', rep('center', ncol(COHORT))), font = list(color = 'white', size = 12)),cells = list(values = rbind(rownames(COHORT),t(as.matrix(unname(COHORT)))), line = list(color = '#506784'),fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222,222,249, 0.65)')),font = list(family = "Arial", size = 12, color = c('#506784'))))
plot_ly(type = 'table', header = list(values = c("<b>Samples</b>","<b>Admin","</b><b>Hannover</b>","<b>Munchen</b>","<b>NULL</b>","<b>Tubingen</b>",names(SAMPLE)),line = list(color = '#506784'),  fill = list(color = 'rgb(235, 100,100,100, 230)'), align = c('Center', rep('center', ncol(SAMPLE))), font = list(color = 'white', size = 12)),cells = list(values = rbind(rownames(SAMPLE),t(as.matrix(unname(SAMPLE)))), line = list(color = '#506784'),fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222,222, 249, 0.65)')),font = list(family = "Arial", size = 12, color = c('#506784'))))
plot_ly(type = 'table', header = list(values = c("<b>Samples</b>","<b>MAL","</b><b>TBC</b>","<b>IICH</b>","<b>NULL</b>",names(STUDY)),line = list(color = '#506784'),  fill = list(color = 'rgb(235, 100,100, 230)'), align = c('Center', rep('center', ncol(STUDY))), font = list(color = 'white', size = 12)),cells = list(values = rbind(rownames(STUDY),t(as.matrix(unname(STUDY)))), line = list(color = '#506784'),fill = list(color = c('rgb(235, 193, 238)', 'rgba(228, 222, 249, 0.65)')),font = list(family = "Arial", size = 12, color = c('#506784'))))

# Unique columns 

table(unique(g$PATIENTID))

occurance <- data.frame(table(g$PATIENTID))
occurance[occurance$Freq > 1,]

# Date format change in various column

d$SAMPLINGDATE <- as.Date(d$SAMPLINGDATE)
f$SAMPLINGDATE <- as.Date(f$SAMPLINGDATE)
f$`Processing Time` <- as.Date(f$`Processing  Time`)
g$`Sampling Time` <- as.Date(g$`Sampling`)
g$`Processing Time` <- as.Date(g$`Processing`)
g$`Secondary DIAGDATE` <- as.Date(g$`Secondary DIAGDATE`)
g$DIAGDATE <- as.Date(g$DIAGDATE)
d$REPOSITIONDATE <- as.Date(d$REPOSITIONDATE)

f$REPOSITIONDATE <- as.Date(as_datetime(f$REPOSITIONDATE))
class(f$REPOSITIONDATE)

am <- c(d$SAMPLINGDATE,f$SAMPLINGDATE,g$`Sampling`,g$`Processing`)
SAM <- as.Date(am, format = "%Y-%m-%d %h:%m")

TIME <- c(d$SAMPLINGDATE,f$SAMPLINGDATE,f$REPOSITIONDATE)
TIME <- as.Date(TIME, format = "%Y-%m-%d-%h-%m")

# Time to Centrifuge calculation , TTC = Time to centrifuge

TTC <- data.frame(SAMPLINGDATE= f$SAMPLINGDATE, PROCESSING= f$`Processing  Time`, Samples = f$SAMPLETYPE , Sites =f$SAMPLELOCATION)
TTC$hours 
head(TTC)
TTC$PROCESSING <- as.Date(as_datetime(TTC$PROCESSING))
TTC$hours <- with(TTC, (difftime(SAMPLINGDATE,PROCESSING,units="hours") ))
TTC

f$`Centrifugation duration`
head(f)
f$PROCESSING <- as.Date(as_datetime(f$`Processing  Time`))
f$`Centrifugation duration` <- with(f, abs(difftime(SAMPLINGDATE,PROCESSING ,units="hours") ))
f

#Time to freeze calculation TTF = Time to freeze

TTF<- data.frame(SAMPLING= f$SAMPLINGDATE,  REPOSITION= f$REPOSITIONDATE , Samples = f$SAMPLETYPE , Sites =f$SAMPLELOCATION)
TTF$hours
head(TTF)

TTF$hours <- with(TTF, abs(difftime(SAMPLING,REPOSITION,units="hours")) )
TTF

f$`Freezing duration`
head(f)
f$PROCESSING <- as.Date(as_datetime(f$`Processing  Time`))
f$`Freezing duration` <- with(f, abs(difftime(SAMPLINGDATE, REPOSITIONDATE,units="hours") ))
f

#Time To centrifuge(TTC) Boxplot 

TTC_new <- TTC[TTC$hours<=12,]
ggplot(TTC_new, aes(x=hours, y=Sites))+ geom_boxplot(outlier.colour="blue", outlier.shape=1,outlier.size=8, notch=TRUE ) +  ggtitle("Centrifugation average duration at Partner sites") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(TTC_new, aes(x=hours, y=Samples))+ geom_boxplot(outlier.colour="blue", outlier.shape=1,outlier.size=8, notch=TRUE ) +  ggtitle("Centrifugation average duration of differernt samples types") +  theme(plot.title = element_text(hjust = 0.5))
  
#Time To Freeze (TTF) Boxplot:

TTF_new <- TTF[TTF$hours<=12,]
ggplot(TTF_new, aes(x=hours, y=Sites))+ geom_boxplot(outlier.colour="red", outlier.shape=1,outlier.size=4, notch=TRUE ) +  ggtitle("Freezing average duration at Partner sites") +  theme(plot.title = element_text(hjust = 0.5))
ggplot(TTF_new, aes(x=hours, y=Samples))+ geom_boxplot(outlier.colour="red", outlier.shape=1,outlier.size=8, notch=TRUE ) +  ggtitle("Samples Freezing average duration of various samples types") +  theme(plot.title = element_text(hjust = 0.5))

eval(parse(dashboard.r, encoding="UTF-8"))
output$plot10 <- renderPlot({
   LMUcountDF %>% ggplot(mapping = aes(x = Standard, y = percentage)) +geom_tile(mapping = aes(fill = Standard))+geom_text(aes(label = `Null count`)) +ggtitle("Null Count of M?nchen (LMU) TX") + theme_light() +theme(plot.title = element_text(hjust = 0.5)) + scale_fill_viridis_d(option = "D")}) 

# Create CSV for export 

write.csv(d,"C:\\Users\\Bhatt\\Desktop\\M2_internsip\\generated files\\New Aliquots.csv", row.names = FALSE)
write.csv(e,"C:\\Users\\Bhatt\\Desktop\\M2_internsip\\generated files\\New Patients.csv", row.names = FALSE)
write.csv(f,"C:\\Users\\Bhatt\\Desktop\\M2_internsip\\generated files\\New samples.csv", row.names = FALSE)
write.csv(g,"C:\\Users\\Bhatt\\Desktop\\M2_internsip\\generated files\\New visits.csv", row.names = FALSE)