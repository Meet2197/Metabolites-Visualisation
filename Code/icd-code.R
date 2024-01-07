
# ICD code extract

NASH<-subset(hesin_diag, diag_icd10=="K758")
MAFLD<-subset(hesin_diag, diag_icd10=="K760")

# DUplication removal for ICD codes

NASH<- as.data.frame(NASH[!duplicated(NASH$eid),])
MAFLD<- as.data.frame(MAFLD[!duplicated(MAFLD$eid),])

# changing value to 1 for designated ICD code disease name : 
MAFLD$MAFLD <-1
NASH$NASH <-1

# icd code split
MAFLD<-select(MAFLD, c(eid_1, MAFLD))

#NASH

NASH<-select(NASH, c(eid_1,NASH))


setnames(NASH, old = colnames(NASH), new = c('eid', 'NASH'))
setnames(MAFLD, old = colnames(MAFLD), new = c('eid', 'MAFLD'))

baseline_nash<-merge(baseline_df,NASH,by="eid",all.x=TRUE, all.y = TRUE)
baseline_nash$NASH[is.na(baseline_nash$NASH)] <- 0

baseline_MAFLD<-merge(baseline_df,MAFLD,by="eid",all.x=TRUE, all.y = TRUE)
baseline_MAFLD$MAFLD[is.na(baseline_MAFLD$MAFLD)] <- 0

#Death Info : 
# This analysis provides information about liver associated death in Liver patients. 

#Extract only K70
K70<-subset(hesin_diag,  startsWith(as.character(hesin_diag$diag_icd10), 'K70'))

#remove duplicates
K70<- as.data.frame(K70[!duplicated(K70$eid),])

# create a new column that K70 equals 1
K70$K70<-1

#select only the ID and K70
K70<-select(K70, c(eid_1, K70))



